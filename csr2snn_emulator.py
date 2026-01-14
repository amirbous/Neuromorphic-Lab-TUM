import scipy as sp
import scipy.sparse
import scipy.sparse.linalg
import numpy as np
import sys
import os
import matplotlib.pyplot as plt
import brian2 as b2

# Ensure we use standard float64 for accuracy matching the numpy version
b2.prefs.core.default_float_dtype = np.float64

# ------------------------------------------------------------------
# 1. HELPER: CONVERT SPARSE MATRIX TO BRIAN2 SYNAPSES
# ------------------------------------------------------------------
def create_synapses_from_sparse(source_group, target_group, sparse_matrix, on_pre_code, variable_name='w'):
    """
    Efficiently creates Brian2 synapses from a scipy.sparse matrix.
    Explicitly removes zero-weight connections to keep the network sparse.
    """
    # 1. Remove explicit zeros stored in the sparse structure
    sparse_matrix.eliminate_zeros()
    
    # 2. (Optional but recommended) Prune numerical noise
    #    Remove weights effectively 0 (e.g. < 1e-12) that prevent sparsity
    sparse_matrix.data[np.abs(sparse_matrix.data) < 1e-12] = 0
    sparse_matrix.eliminate_zeros()

    # 3. Convert to COO for iteration
    coo = sparse_matrix.tocoo()
    
    # 4. Define Synapses object
    S = b2.Synapses(source_group, target_group, model=f'{variable_name} : 1', on_pre=on_pre_code)
    
    # 5. Connect efficiently using array passing
    S.connect(i=coo.row, j=coo.col)
    
    # 6. Assign weights
    setattr(S, variable_name, coo.data)
    
    return S

def run_neurofem_brian2(GTAG, omega_f, threshs, v_reset, input_bias,
                        gamma, dt_val, n_timesteps,
                        lambda_d, lambda_v, kp, ki, sigma_v):
    """
    Replaces the manual numpy loop with a Brian2 network simulation.
    """
    
    # --- 1. Setup Simulation Environment ---
    b2.start_scope()
    b2.defaultclock.dt = dt_val * b2.second

    n_neurons = GTAG.shape[0]
    n_sys = gamma.shape[0]

    # --- 2. Define Main Spiking Population ---
    
    # Noise Scaling for Euler integration
    sigma_sde = sigma_v / np.sqrt(dt_val)

    # Differential equations
    eqs_neurons = '''
    dv/dt = (-lambda_v * v + kp * (u1 + bias) + ki * u_int + u2) / second + sigma_sde * xi * second**-0.5 : 1
    du1/dt = -lambda_d * u1 / second : 1
    du2/dt = -lambda_d * u2 / second : 1
    du_int/dt = (u1 + bias) / second : 1
    
    bias : 1
    v_th : 1
    v_rst : 1
    '''

    # Create the NeuronGroup
    neurons = b2.NeuronGroup(n_neurons, eqs_neurons, 
                             threshold='v > v_th', 
                             reset='v = v_rst', 
                             method='euler',
                             namespace={'kp': kp, 'ki': ki, 'lambda_v': lambda_v, 
                                        'lambda_d': lambda_d, 'sigma_sde': sigma_sde})
    
    # --- Initialize State Variables ---
    neurons.v = 0
    neurons.u1 = 0 
    neurons.u2 = 0
    neurons.u_int = 0
    neurons.v_th = threshs
    neurons.v_rst = v_reset
    neurons.bias = input_bias

    # --- 3. Define Connectivity (Recurrent) ---
    # The helper now strips zeros, resulting in drastically fewer synapses
    S_GTAG = create_synapses_from_sparse(neurons, neurons, GTAG, 
                                         on_pre_code='u1 += w', variable_name='w')

    S_omega_u2 = create_synapses_from_sparse(neurons, neurons, omega_f, 
                                             on_pre_code='u2 += lambda_d * w', variable_name='w')

    S_omega_v = create_synapses_from_sparse(neurons, neurons, omega_f, 
                                            on_pre_code='v -= w', variable_name='w')

    # --- 4. Define Readout Population ---
    eqs_readout = '''
    dx/dt = -lambda_d * x / second : 1
    '''
    readout = b2.NeuronGroup(n_sys, eqs_readout, method='euler', 
                             namespace={'lambda_d': lambda_d})
    readout.x = 0
    
    S_readout = create_synapses_from_sparse(neurons, readout, gamma.T, 
                                            on_pre_code='x += w', variable_name='w')

    # --- 5. Monitors ---
    mon_readout = b2.StateMonitor(readout, 'x', record=True, dt=dt_val*b2.second)
    mon_spikes = b2.SpikeMonitor(neurons)

    # --- 6. Run ---
    total_synapses = S_GTAG.N + S_omega_u2.N + S_readout.N + S_omega_v.N
    print(f"  [Brian2] Building network (N={n_neurons}, Synapses={total_synapses} (optimized))...")
    
    b2.run(n_timesteps * dt_val * b2.second)
    
    return mon_readout.x, mon_spikes

# ------------------------------------------------------------------
# 3. EXISTING LOADING & MATRIX GEN (Unchanged)
# ------------------------------------------------------------------
def load_csr_matrix(model_name):
    mtx_file = f"data/matrix/{model_name}_mtx.txt"
    rhs_file = f"data/matrix/{model_name}_rhs.txt"
    
    # Load raw data
    raw = np.loadtxt(mtx_file, skiprows=1)
    with open(mtx_file, 'r') as f:
        head = f.readline().split()
    n_rows, n_cols = int(head[0]), int(head[1])
    
    rows = raw[:, 0].astype(int) - 1
    cols = raw[:, 1].astype(int) - 1
    vals = raw[:, 2]
    
    A = sp.sparse.csr_matrix((vals, (rows, cols)), shape=(n_rows, n_cols))
    b_vec = np.loadtxt(rhs_file, skiprows=1)
    return A, b_vec

def generate_gamma_sparse(n_sys, neurons_per_mesh_point, gamma_norm):
    n_neurons = neurons_per_mesh_point * n_sys
    data, rows, cols = [], [], []
    for pt in range(n_sys):
        start = pt * neurons_per_mesh_point
        half = neurons_per_mesh_point // 2
        rows.extend(neurons_per_mesh_point * [pt])
        cols.extend(range(start, start + neurons_per_mesh_point))
        data.extend(half * [-gamma_norm] + half * [gamma_norm])
    return sp.sparse.csr_array((np.array(data), (np.array(rows), np.array(cols))),
                               shape=(n_sys, n_neurons))

def create_spiking_fem_network_sparse(A_sys, gamma, lambda_d, lambda_v,
                                      mu, nu, tau_A, gamma_norm):
    n_sys, n_neurons = gamma.shape
    gTg = gamma.T @ gamma
    
    # Identity matrix
    I = sp.sparse.eye_array(n_neurons)
    
    # slow weight matrix
    omega_f = gTg + mu * (lambda_d**2) * I
    
    # fast weight matrix
    GTAG = (gamma.T @ (A_sys / tau_A) @ gamma).tocsc()
    
    # thresholds and reset
    threshs = 0.5 * (nu * lambda_d**2 + (gamma_norm**2) * np.ones(n_neurons))
    v_reset = threshs - (gamma_norm**2 + mu * (lambda_d**2)) * np.ones(n_neurons)
    
    # subtract diagonal self-reset term from omega_f
    omega_f = omega_f - (gamma_norm**2 + mu * (lambda_d**2)) * I
    
    # IMPORTANT: Clean matrices here too
    omega_f.eliminate_zeros()
    GTAG.eliminate_zeros()
    
    return omega_f, GTAG, threshs, v_reset

# ------------------------------------------------------------------
# 4. MAIN SOLVER
# ------------------------------------------------------------------
def neurofem_solve_brian(A, b_vec, neurons_per_mesh_point=16, sim_time=2.0):
    n = A.shape[0]

    # Precondition
    diag = A.diagonal()
    diag[diag == 0] = 1.0
    Pinv = sp.sparse.diags(1.0 / diag)
    A_sys = Pinv @ A
    b_sys = Pinv @ b_vec

    # Parameters
    gamma_pow2 = -6
    s_gamma = 7 - gamma_pow2
    gamma_norm = (2**gamma_pow2 - 2**-s_gamma)
    omega_s_pow2 = 1
    omega_max = 2**omega_s_pow2 - 2**-(7 - omega_s_pow2)
    lambda_d = 8.0
    lambda_v = 16.0
    ki = 16.0
    kp = 4.0
    dt = 2**(-12)
    sigma_v = 0.00225
    mu = 0.0
    nu = 0.0

    n_timesteps = int(sim_time / dt)
    print(f"Timesteps: {n_timesteps}, dt: {dt}")

    # Scaling with safety factor for saturation
    A_max = np.amax(np.abs(A_sys))
    tau_A = ((gamma_norm**2) * A_max / omega_max) * 4.0 

    # Network Gen
    gamma = generate_gamma_sparse(n, neurons_per_mesh_point, gamma_norm)
    omega_f, GTAG, threshs, v_reset = create_spiking_fem_network_sparse(
        -A_sys, gamma, lambda_d, lambda_v, mu, nu, tau_A, gamma_norm
    )

    bias = gamma.T @ (b_sys / tau_A)

    print(f"  [Parameters] dt = {dt:.5f}s, Total time = {sim_time}s")
    
    # --- RUN BRIAN2 SIMULATION ---
    output, spike_mon = run_neurofem_brian2(
        GTAG, omega_f, threshs, v_reset, bias,
        gamma, dt, n_timesteps,
        lambda_d, lambda_v, kp, ki, sigma_v
    )

    # Decode
    burn_in = int(0.2 * output.shape[1])
    x_snn = np.mean(output[:, burn_in:], axis=1)
    return x_snn

# ------------------------------------------------------------------
# 5. MAIN EXECUTION
# ------------------------------------------------------------------
def main():
    model_name = sys.argv[1] if len(sys.argv) > 1 else "Sphere_00"
    print(f"=== NeuroFEM w/ Brian2 Emulator ===")
    print(f"Problem: {model_name}")

    A, b_vec = load_csr_matrix(model_name)
    
    # Solve
    x_snn = neurofem_solve_brian(A, b_vec, neurons_per_mesh_point=16, sim_time=2.0)

    # Compare
    x_exact = sp.sparse.linalg.spsolve(A, b_vec)
    
    # Calculate metrics
    corr = np.corrcoef(x_exact, x_snn)[0, 1]
    mae = np.mean(np.abs(x_exact - x_snn))
    
    print(f"\n[Results]")
    print(f"  Correlation: {corr:.4f}")
    print(f"  MAE:         {mae:.4f}")

    # Plot
    plt.figure(figsize=(10, 5))
    plt.subplot(121)
    plt.plot(x_exact, label='Exact')
    plt.plot(x_snn, label='Brian2 SNN', alpha=0.7)
    plt.legend()
    plt.title('Solution Profile')
    
    plt.subplot(122)
    plt.scatter(x_exact, x_snn, alpha=0.5)
    plt.plot([x_exact.min(), x_exact.max()], [x_exact.min(), x_exact.max()], 'k--')
    plt.title(f'Correlation (r={corr:.3f})')
    
    plt.tight_layout()
    plt.savefig(f"{model_name}_neurofem_brian2_results.png", dpi=300)

if __name__ == "__main__":
    main()