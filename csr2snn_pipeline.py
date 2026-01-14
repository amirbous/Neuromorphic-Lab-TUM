import scipy as sp
import scipy.sparse
import scipy.sparse.linalg
import numpy as np
import sys
import os
import matplotlib.pyplot as plt

# SpiNNaker2 imports
try:
    from spinnaker2 import snn, hardware
    from spinnaker2.experiment_backends import BackendSettings
    from spinnaker2.experiment_backends.backend_settings import ROUTING
    HAS_SPINNAKER2 = True
except ImportError:
    HAS_SPINNAKER2 = False
    print("WARNING: SpiNNaker2 not available")

# For validation only (Brian2 emulator)
try:
    import brian2 as b2
    HAS_BRIAN2 = True
except ImportError:
    HAS_BRIAN2 = False
    print("WARNING: Brian2 not available for validation")

# ==========================================
# 1. EXACT SAME LOADING & NEUROFEM FUNCTIONS
# ==========================================
def load_csr_matrix(model_name):
    mtx_file = f"data/matrix/{model_name}_mtx.txt"
    rhs_file = f"data/matrix/{model_name}_rhs.txt"
    
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
    
    omega_f.eliminate_zeros()
    GTAG.eliminate_zeros()
    
    return omega_f, GTAG, threshs, v_reset

# ==========================================
# 2. SPINNAKER2 HARDWARE CONVERTER
# ==========================================
def convert_to_spinnaker2(GTAG, omega_f, threshs, v_reset, bias, gamma):
    """
    Convert NeuroFEM floating-point matrices to SpiNNaker2 integer format.
    This is the ONLY adaptation needed for hardware.
    """
    n_neurons = GTAG.shape[0]
    
    # Hardware constraints
    MAX_SYNAPSE_WEIGHT = 15      # 4-bit signed: -15 to 15
    MAX_NEURON_PARAM = 32767     # 15-bit signed
    HARDWARE_THRESHOLD = 1000    # Fixed threshold for all neurons
    
    # 1. Combine weight matrices for hardware (single weight matrix)
    # We'll use GTAG as primary, omega_f scaled down
    alpha = 0.1  # Balance factor
    W_combined = GTAG + alpha * omega_f
    
    # 2. Scale weights to hardware range (-15 to 15)
    max_weight = np.max(np.abs(W_combined.data))
    if max_weight > 0:
        weight_scale = 0.7 * MAX_SYNAPSE_WEIGHT / max_weight
    else:
        weight_scale = 1.0
    
    W_scaled = W_combined * weight_scale
    
    # Convert to integer and clip
    W_coo = W_scaled.tocoo()
    W_coo.data = np.round(W_coo.data).astype(np.int32)
    W_coo.data = np.clip(W_coo.data, -MAX_SYNAPSE_WEIGHT, MAX_SYNAPSE_WEIGHT)
    
    # Build connections list for SpiNNaker2
    connections = []
    for i, j, w in zip(W_coo.row, W_coo.col, W_coo.data):
        if abs(w) > 0:  # Skip zero weights
            connections.append([int(i), int(j), int(w), 1])  # [source, target, weight, delay]
    
    # 3. Scale biases to hardware range
    bias_scaled = bias * weight_scale
    bias_scaled = np.round(bias_scaled).astype(np.int32)
    bias_scaled = np.clip(bias_scaled, -MAX_NEURON_PARAM, MAX_NEURON_PARAM)
    
    # 4. Hardware uses fixed threshold
    threshold_hw = HARDWARE_THRESHOLD
    
    print(f"  [SpiNNaker2] Converted {n_neurons} neurons, {len(connections)} synapses")
    print(f"  [SpiNNaker2] Weight range: [{W_coo.data.min()}, {W_coo.data.max()}]")
    print(f"  [SpiNNaker2] Bias range: [{bias_scaled.min()}, {bias_scaled.max()}]")
    print(f"  [SpiNNaker2] Fixed threshold: {threshold_hw}")
    
    return connections, bias_scaled, threshold_hw, weight_scale

# ==========================================
# 3. SPINNAKER2 HARDWARE RUNNER
# ==========================================
def run_on_spinnaker2(connections, bias, threshold, n_neurons, timesteps=5000):
    """Run network on SpiNNaker2 hardware"""
    if not HAS_SPINNAKER2:
        raise RuntimeError("SpiNNaker2 not available")
    
    print(f"  [Hardware] Deploying to SpiNNaker2...")
    
    # Build SpiNNaker2 network
    lif_params = {
        "threshold": threshold,
        "i_offset": bias.astype(np.int32),
        "reset": "reset_by_subtraction",
        "alpha_decay": 1.0,  # No leak
        "beta_decay": 0.0,
        "v_reset": 0
    }
    
    pop = snn.Population(n_neurons, "lif", lif_params,
                        name="neurofem_pop", record=["spikes"])
    pop.set_max_atoms_per_core(32)  # Distribute across cores
    
    proj = snn.Projection(pop, pop, connections)
    net = snn.Network("NeuroFEM_HW")
    net.add(pop, proj)
    
    # Configure hardware
    cfg = BackendSettings()
    cfg.switch_serial = False
    cfg.routing_type = ROUTING.C2C
    cfg.rebuild_apps = True
    
    # Connect to hardware
    try:
        # Try multi-chip board first
        ip = os.environ.get("STM_IP", "192.168.1.2")
        print(f"    Connecting to SpiNNaker2 at {ip}...")
        hw = hardware.SpiNNcloud48NodeBoard(stm_ip=ip)
    except Exception as e:
        # Fall back to single chip
        ip = os.environ.get("S2_IP", "192.168.1.17")
        print(f"    Connecting to single chip at {ip}...")
        hw = hardware.SpiNNaker2Chip(eth_ip=ip)
    
    # Run on hardware
    hw.run(net, timesteps, experiment_backend_settings=cfg)
    
    # Get spikes
    spikes = pop.get_spikes()
    
    print(f"  [Hardware] Total spikes: {sum(len(s) for s in spikes.values())}")
    
    return spikes
# ==========================================
# 4. BRIAN2 EMULATOR (FIXED)
# ==========================================
def run_brian2_emulator(GTAG, omega_f, threshs, v_reset, input_bias,  # <--- RENAMED HERE
                        gamma, dt, n_timesteps,
                        lambda_d, lambda_v, kp, ki, sigma_v):
    """Run NeuroFEM on Brian2 for validation"""
    if not HAS_BRIAN2:
        raise RuntimeError("Brian2 not available")
    
    b2.start_scope()
    b2.defaultclock.dt = dt * b2.second
    
    n_neurons = GTAG.shape[0]
    n_sys = gamma.shape[0]
    
    # Brian2 implementation (same as original)
    sigma_sde = sigma_v / np.sqrt(dt)
    
    # 'bias' inside the string refers to the NeuronGroup state variable
    eqs_neurons = '''
    dv/dt = (-lambda_v * v + kp * (u1 + bias) + ki * u_int + u2) / second + sigma_sde * xi * second**-0.5 : 1
    du1/dt = -lambda_d * u1 / second : 1
    du2/dt = -lambda_d * u2 / second : 1
    du_int/dt = (u1 + bias) / second : 1
    
    bias : 1
    v_th : 1
    v_rst : 1
    '''
    
    neurons = b2.NeuronGroup(n_neurons, eqs_neurons, 
                            threshold='v > v_th', 
                            reset='v = v_rst', 
                            method='euler',
                            namespace={'kp': kp, 'ki': ki, 'lambda_v': lambda_v, 
                                      'lambda_d': lambda_d, 'sigma_sde': sigma_sde})
    
    neurons.v = 0
    neurons.u1 = 0 
    neurons.u2 = 0
    neurons.u_int = 0
    neurons.v_th = threshs
    neurons.v_rst = v_reset
    
    # Assign the Python argument (input_bias) to the Brian2 state variable (bias)
    neurons.bias[:] = input_bias   # <--- UPDATED ASSIGNMENT
    
    # Helper function to create synapses
    def create_synapses(source, target, matrix, on_pre_code, var_name='w'):
        matrix.eliminate_zeros()
        coo = matrix.tocoo()
        S = b2.Synapses(source, target, model=f'{var_name} : 1', on_pre=on_pre_code)
        S.connect(i=coo.row, j=coo.col)
        setattr(S, var_name, coo.data)
        return S
    
    # Create connections
    S_GTAG = create_synapses(neurons, neurons, GTAG, 'u1 += w')
    S_omega_u2 = create_synapses(neurons, neurons, omega_f, 'u2 += lambda_d * w')
    S_omega_v = create_synapses(neurons, neurons, omega_f, 'v -= w')
    
    # Readout
    eqs_readout = 'dx/dt = -lambda_d * x / second : 1'
    readout = b2.NeuronGroup(n_sys, eqs_readout, method='euler', 
                            namespace={'lambda_d': lambda_d})
    readout.x = 0
    S_readout = create_synapses(neurons, readout, gamma.T, 'x += w')
    
    # Monitors
    mon_readout = b2.StateMonitor(readout, 'x', record=True, dt=dt*b2.second)
    mon_spikes = b2.SpikeMonitor(neurons)
    
    print(f"  [Brian2] Running {n_timesteps} timesteps...")
    b2.run(n_timesteps * dt * b2.second)
    
    return mon_readout.x, mon_spikes
# ==========================================
# 5. DECODE SOLUTION (SAME FOR BOTH)
# ==========================================
def decode_solution(output, gamma, neurons_per_mesh_point, sim_time, burn_in_frac=0.2):
    """Decode solution from spike counts (same for Brian2 and SpiNNaker2)"""
    n_sys = gamma.shape[0]
    
    if isinstance(output, np.ndarray):  # Brian2 output (state monitor)
        # output is [n_sys x n_timesteps]
        burn_in = int(output.shape[1] * burn_in_frac)
        x_snn = np.mean(output[:, burn_in:], axis=1)
        
    elif isinstance(output, dict):  # SpiNNaker2 spikes
        # output is spike dict {neuron_id: [spike_times]}
        n_neurons = len(output)
        n_timesteps = int(sim_time * 1000)  # Convert to ms
        burn_in = int(n_timesteps * burn_in_frac)
        
        # Count spikes in analysis window
        spike_counts = np.zeros(n_neurons)
        for neuron_id, spike_times in output.items():
            valid_spikes = [t for t in spike_times if t >= burn_in and t < n_timesteps]
            spike_counts[neuron_id] = len(valid_spikes)
        
        # Decode using gamma matrix
        x_snn = np.zeros(n_sys)
        gamma_norm = abs(gamma[0, 0])  # Get gamma value
        
        for var_idx in range(n_sys):
            start = var_idx * neurons_per_mesh_point
            half = neurons_per_mesh_point // 2
            
            # Count spikes for this variable
            pos_spikes = 0
            neg_spikes = 0
            
            for i in range(half):
                neuron_id = start + i
                neg_spikes += spike_counts[neuron_id]
            
            for i in range(half, neurons_per_mesh_point):
                neuron_id = start + i
                pos_spikes += spike_counts[neuron_id]
            
            # Decode
            rate = (pos_spikes - neg_spikes) / (n_timesteps - burn_in) * 1000  # Hz
            x_snn[var_idx] = rate * gamma_norm
    
    else:
        raise ValueError(f"Unknown output type: {type(output)}")
    
    return x_snn

# ==========================================
# 6. MAIN SOLVER
# ==========================================
def solve_neurofem(model_name, use_hardware=True, validate=True, 
                   neurons_per_mesh_point=8, sim_time=2.0):
    """
    Solve FEM problem using NeuroFEM.
    
    Parameters:
    -----------
    model_name : str
        Name of the model to load
    use_hardware : bool
        Whether to run on SpiNNaker2 hardware
    validate : bool
        Whether to validate with Brian2 emulator
    neurons_per_mesh_point : int
        Number of neurons per variable (must be even)
    sim_time : float
        Simulation time in seconds
    """
    print(f"\n{'='*70}")
    print(f"NEUROFEM SOLVER: {model_name}")
    print(f"{'='*70}")
    
    # 1. Load the problem
    print("\n[1] Loading linear system...")
    A, b_vec = load_csr_matrix(model_name)
    n = A.shape[0]
    print(f"    Size: {n}×{n}, Nonzeros: {A.nnz}")
    
    # Get exact solution for comparison
    exact_solution = sp.sparse.linalg.spsolve(A, b_vec)
    
    # 2. NeuroFEM parameters (from paper)
    print("\n[2] Setting up NeuroFEM parameters...")
    
    # NeuroFEM parameters (exact same as your Brian2 code)
    gamma_pow2 = -6
    s_gamma = 7 - gamma_pow2
    gamma_norm = (2**gamma_pow2 - 2**-s_gamma)
    
    omega_s_pow2 = 1
    omega_max = 2**omega_s_pow2 - 2**-(7 - omega_s_pow2)
    
    lambda_d = 8.0
    lambda_v = 16.0
    ki = 16.0
    kp = 4.0
    dt = 2**(-12)      # ~0.244 ms
    sigma_v = 0.00225
    mu = 0.0
    nu = 0.0
    
    n_timesteps = int(sim_time / dt)
    print(f"    γ = {gamma_norm:.6f}, dt = {dt:.6f}s, Steps = {n_timesteps}")
    
    # 3. Create NeuroFEM network (exact same as your Brian2 code)
    print("\n[3] Creating NeuroFEM network...")
    
    # No Jacobi preconditioning - using original A and b
    A_sys = A.copy()
    b_sys = b_vec.copy()
    
    # Scaling factor
    A_max = np.amax(np.abs(A_sys))
    tau_A = ((gamma_norm**2) * A_max / omega_max) * 4.0
    
    # Generate gamma matrix and weight matrices
    gamma = generate_gamma_sparse(n, neurons_per_mesh_point, gamma_norm)
    omega_f, GTAG, threshs, v_reset = create_spiking_fem_network_sparse(
        -A_sys, gamma, lambda_d, lambda_v, mu, nu, tau_A, gamma_norm
    )
    
    bias = gamma.T @ (b_sys / tau_A)
    
    print(f"    Neurons: {gamma.shape[1]}, Synapses: {GTAG.nnz + omega_f.nnz}")
    
    # 4. Run Brian2 validation (if requested)
    brian2_solution = None
    if validate and HAS_BRIAN2:
        print("\n[4] Running Brian2 emulator (validation)...")
        try:
            b2_output, b2_spikes = run_brian2_emulator(
                GTAG, omega_f, threshs, v_reset, bias,
                gamma, dt, n_timesteps,
                lambda_d, lambda_v, kp, ki, sigma_v
            )
            
            brian2_solution = decode_solution(
                b2_output, gamma, neurons_per_mesh_point, sim_time
            )
            
            b2_error = exact_solution - brian2_solution
            b2_rel_err = np.linalg.norm(b2_error) / np.linalg.norm(exact_solution) * 100
            b2_corr = np.corrcoef(exact_solution, brian2_solution)[0, 1]
            
            print(f"    Brian2 - Relative error: {b2_rel_err:.2f}%")
            print(f"    Brian2 - Correlation: {b2_corr:.4f}")
            
        except Exception as e:
            print(f"    Brian2 error: {e}")
    
    # 5. Run on SpiNNaker2 hardware (if requested)
    spinnaker2_solution = None
    if use_hardware and HAS_SPINNAKER2:
        print("\n[5] Running on SpiNNaker2 hardware...")
        try:
            # Convert to hardware format
            connections, bias_hw, threshold_hw, weight_scale = convert_to_spinnaker2(
                GTAG, omega_f, threshs, v_reset, bias, gamma
            )
            
            # Run on hardware
            spikes_hw = run_on_spinnaker2(
                connections, bias_hw, threshold_hw,
                n_neurons=gamma.shape[1],
                timesteps=int(sim_time * 1000)  # Convert to ms
            )
            
            # Decode solution
            spinnaker2_solution = decode_solution(
                spikes_hw, gamma, neurons_per_mesh_point, sim_time
            )
            
            # Apply inverse weight scaling
            spinnaker2_solution = spinnaker2_solution / weight_scale
            
            hw_error = exact_solution - spinnaker2_solution
            hw_rel_err = np.linalg.norm(hw_error) / np.linalg.norm(exact_solution) * 100
            hw_corr = np.corrcoef(exact_solution, spinnaker2_solution)[0, 1]
            
            print(f"    SpiNNaker2 - Relative error: {hw_rel_err:.2f}%")
            print(f"    SpiNNaker2 - Correlation: {hw_corr:.4f}")
            
            # Save solution
            os.makedirs("data/sol", exist_ok=True)
            output_file = f"data/sol/{model_name}_spinnaker2.txt"
            np.savetxt(output_file, spinnaker2_solution, fmt="%.8f")
            print(f"    Solution saved to: {output_file}")
            
        except Exception as e:
            print(f"    SpiNNaker2 error: {e}")
    
    # 6. Compare and plot
    print("\n[6] Generating comparison plots...")
    try:
        fig, axes = plt.subplots(2, 2, figsize=(14, 10))
        
        solutions_to_plot = []
        labels = []
        
        # Exact solution
        solutions_to_plot.append(exact_solution)
        labels.append('Exact')
        
        # Brian2 solution
        if brian2_solution is not None:
            solutions_to_plot.append(brian2_solution)
            labels.append('Brian2')
        
        # SpiNNaker2 solution
        if spinnaker2_solution is not None:
            solutions_to_plot.append(spinnaker2_solution)
            labels.append('SpiNNaker2')
        
        # Plot 1: Solution comparison
        ax = axes[0, 0]
        n_plot = min(50, n)
        idx = np.arange(n_plot)
        
        for sol, label in zip(solutions_to_plot, labels):
            ax.plot(idx, sol[:n_plot], '-', label=label, alpha=0.7, 
                   linewidth=2 if label == 'Exact' else 1.5)
        
        ax.set_xlabel('Node Index')
        ax.set_ylabel('Solution Value')
        ax.legend()
        ax.grid(True, alpha=0.3)
        ax.set_title(f'Solution Comparison ({model_name})')
        
        # Plot 2: Correlation plot
        ax = axes[0, 1]
        if brian2_solution is not None:
            ax.scatter(exact_solution, brian2_solution, alpha=0.6, s=20, 
                      label=f'Brian2 (r={np.corrcoef(exact_solution, brian2_solution)[0,1]:.3f})')
        
        if spinnaker2_solution is not None:
            ax.scatter(exact_solution, spinnaker2_solution, alpha=0.6, s=20, marker='s',
                      label=f'SpiNNaker2 (r={np.corrcoef(exact_solution, spinnaker2_solution)[0,1]:.3f})')
        
        lims = [exact_solution.min(), exact_solution.max()]
        ax.plot(lims, lims, 'k--', linewidth=1, label='Perfect')
        ax.set_xlabel('Exact Solution')
        ax.set_ylabel('SNN Solution')
        ax.legend()
        ax.grid(True, alpha=0.3)
        ax.set_title('Correlation Plot')
        
        # Plot 3: Error distribution
        ax = axes[1, 0]
        errors = []
        error_labels = []
        
        if brian2_solution is not None:
            errors.append(exact_solution - brian2_solution)
            error_labels.append('Brian2')
        
        if spinnaker2_solution is not None:
            errors.append(exact_solution - spinnaker2_solution)
            error_labels.append('SpiNNaker2')
        
        for err, label in zip(errors, error_labels):
            ax.hist(err, bins=30, alpha=0.5, label=label, edgecolor='black')
        
        ax.axvline(0, color='k', linestyle='--', linewidth=1)
        ax.set_xlabel('Error')
        ax.set_ylabel('Frequency')
        ax.legend()
        ax.grid(True, alpha=0.3)
        ax.set_title('Error Distribution')
        
        # Plot 4: Summary statistics
        ax = axes[1, 1]
        ax.axis('off')
        
        stats_text = "Solution Statistics:\n"
        stats_text += "=" * 40 + "\n"
        stats_text += f"Exact Solution:\n"
        stats_text += f"  Mean: {np.mean(exact_solution):.4f}\n"
        stats_text += f"  Std:  {np.std(exact_solution):.4f}\n"
        stats_text += f"  Range: [{exact_solution.min():.4f}, {exact_solution.max():.4f}]\n\n"
        
        if brian2_solution is not None:
            b2_error = exact_solution - brian2_solution
            stats_text += f"Brian2 Emulator:\n"
            stats_text += f"  Rel. Error: {np.linalg.norm(b2_error)/np.linalg.norm(exact_solution)*100:.2f}%\n"
            stats_text += f"  MAE: {np.mean(np.abs(b2_error)):.4f}\n"
            stats_text += f"  Corr: {np.corrcoef(exact_solution, brian2_solution)[0,1]:.4f}\n\n"
        
        if spinnaker2_solution is not None:
            hw_error = exact_solution - spinnaker2_solution
            stats_text += f"SpiNNaker2 Hardware:\n"
            stats_text += f"  Rel. Error: {np.linalg.norm(hw_error)/np.linalg.norm(exact_solution)*100:.2f}%\n"
            stats_text += f"  MAE: {np.mean(np.abs(hw_error)):.4f}\n"
            stats_text += f"  Corr: {np.corrcoef(exact_solution, spinnaker2_solution)[0,1]:.4f}\n"
        
        ax.text(0.1, 0.5, stats_text, fontfamily='monospace', 
               verticalalignment='center', fontsize=10)
        
        plt.tight_layout()
        plot_file = f"{model_name}_comparison.png"
        plt.savefig(plot_file, dpi=150, bbox_inches='tight')
        print(f"    Plot saved to: {plot_file}")
        plt.close()
        
    except Exception as e:
        print(f"    Plotting error: {e}")
    
    print(f"\n{'='*70}")
    print("SOLUTION COMPLETE!")
    print(f"{'='*70}")
    
    return {
        'exact': exact_solution,
        'brian2': brian2_solution,
        'spinnaker2': spinnaker2_solution
    }

# ==========================================
# 7. COMMAND LINE INTERFACE
# ==========================================
def main():
    """Command line interface"""
    if len(sys.argv) < 2:
        print("Usage: python neurofem_spinnaker2.py <model_name> [options]")
        print("\nOptions:")
        print("  --no-hardware    Run Brian2 emulator only (no hardware)")
        print("  --no-validate    Skip Brian2 validation")
        print("  --neurons N      Neurons per variable (default: 8)")
        print("\nExamples:")
        print("  python neurofem_spinnaker2.py Sphere_00")
        print("  python neurofem_spinnaker2.py Sphere_00 --no-hardware")
        print("  python neurofem_spinnaker2.py Sphere_00 --neurons 16")
        sys.exit(1)
    
    model_name = sys.argv[1]
    
    # Parse options
    use_hardware = '--no-hardware' not in sys.argv
    validate = '--no-validate' not in sys.argv
    
    # Parse neurons per variable
    neurons_per_mesh_point = 16
    for i, arg in enumerate(sys.argv):
        if arg == '--neurons' and i + 1 < len(sys.argv):
            neurons_per_mesh_point = int(sys.argv[i + 1])
    
    # Check hardware availability
    if use_hardware and not HAS_SPINNAKER2:
        print("WARNING: SpiNNaker2 hardware not available. Using Brian2 only.")
        use_hardware = False
    
    # Check Brian2 availability
    if validate and not HAS_BRIAN2:
        print("WARNING: Brian2 not available for validation.")
        validate = False
    
    # Run solver
    results = solve_neurofem(
        model_name=model_name,
        use_hardware=use_hardware,
        validate=validate,
        neurons_per_mesh_point=neurons_per_mesh_point,
        sim_time=2.0
    )
    
    # Print summary
    print("\nFINAL RESULTS:")
    print("-" * 40)
    
    if results['brian2'] is not None:
        b2_err = np.linalg.norm(results['exact'] - results['brian2']) / np.linalg.norm(results['exact']) * 100
        print(f"Brian2 emulator: {b2_err:.2f}% error")
    
    if results['spinnaker2'] is not None:
        hw_err = np.linalg.norm(results['exact'] - results['spinnaker2']) / np.linalg.norm(results['exact']) * 100
        print(f"SpiNNaker2 hardware: {hw_err:.2f}% error")
    
    print("-" * 40)

if __name__ == "__main__":
    main()