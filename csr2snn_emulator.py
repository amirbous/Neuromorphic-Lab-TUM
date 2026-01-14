import scipy as sp
import scipy.sparse
import scipy.sparse.linalg
import numpy as np
import sys
import os
import matplotlib.pyplot as plt

# ------------------------------------------------------------------
# 1. LOAD THE USER’S LINEAR SYSTEM
# ------------------------------------------------------------------
def load_csr_matrix(model_name):
    """Load the sparse matrix A and right‑hand side b from the provided files."""
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
    b = np.loadtxt(rhs_file, skiprows=1)
    return A, b

# ------------------------------------------------------------------
# 2. GENERATE THE Γ MATRIX (exactly as in the notebook)
# ------------------------------------------------------------------
def generate_gamma_sparse(n_sys, neurons_per_mesh_point, gamma_norm):
    """
    Creates the sparse readout kernel Γ.
    For each mesh variable, half the neurons have weight +γ, half have weight –γ.
    """
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

# ------------------------------------------------------------------
# 3. BUILD THE NEUROFEM NETWORK (exactly as in the notebook)
# ------------------------------------------------------------------
def create_spiking_fem_network_sparse(A_sys, gamma, lambda_d, lambda_v,
                                      mu, nu, tau_A, gamma_norm):
    """
    Returns the sparse weight matrices, thresholds, and reset voltages
    as defined in the NeuroFEM notebook.
    """
    n_sys, n_neurons = gamma.shape
    gTg = gamma.T @ gamma
    # slow weight matrix ω_f
    omega_f = gTg + mu * (lambda_d**2) * sp.sparse.eye_array(n_neurons)
    # fast weight matrix G^T A G / τ_A
    GTAG = (gamma.T @ (A_sys / tau_A) @ gamma).tocsc()
    # thresholds and reset voltages
    threshs = 0.5 * (nu * lambda_d**2 + (gamma_norm**2) * np.ones(n_neurons))
    v_reset = threshs - (gamma_norm**2 + mu * (lambda_d**2)) * np.ones(n_neurons)
    omega_f = omega_f - (gamma_norm**2 + mu * (lambda_d**2)) * sp.sparse.eye_array(n_neurons)
    return omega_f, GTAG, threshs, v_reset

# ------------------------------------------------------------------
# 4. SIMULATE THE NEUROFEM DYNAMICS (exactly as in the notebook)
# ------------------------------------------------------------------
def run_neurofem_simulation(GTAG, omega_f, threshs, v_reset, bias,
                            gamma, dt, n_timesteps,
                            lambda_d, lambda_v, kp, ki, sigma_v):
    """
    Simulates the NeuroFEM dynamics using the exact update equations
    from the notebook (function `run_model`).
    """
    n_neurons = GTAG.shape[0]
    n_sys = gamma.shape[0]

    # state variables
    V = np.zeros(n_neurons)
    u1 = np.zeros(n_neurons)
    u2 = np.zeros(n_neurons)
    u_int = np.zeros(n_neurons)
    output = np.zeros((n_sys, n_timesteps))
    spikes = np.zeros((n_neurons, 1), dtype=np.int8)

    for step in range(1, n_timesteps):
        # previous spikes
        spk = spikes[:, 0].astype(float)

        # update slow variables (eqs. 6–7)
        du1 = dt * (-lambda_d * u1) + GTAG @ spk
        u1 += du1
        du2 = dt * (-lambda_d * u2) + lambda_d * omega_f @ spk
        u2 += du2

        # error and integral (eq. 8)
        u_err = u1 + bias
        u_int += dt * u_err

        # membrane potential (eq. 9)
        dv = (dt * (-lambda_v * V + kp * u_err + ki * u_int + u2)
              - omega_f.dot(spk)
              + sigma_v * np.random.randn(n_neurons))
        V += dv

        # spike generation and reset
        spikes[:, 0] = (V > threshs).astype(np.int8)
        V[spikes[:, 0] > 0] = v_reset[spikes[:, 0] > 0]

        # low‑pass filtered readout (eq. 5)
        output[:, step] = (1 - dt * lambda_d) * output[:, step - 1] + gamma.dot(spikes[:, 0])

    return output, spikes

# ------------------------------------------------------------------
# 5. MAIN SOLVER
# ------------------------------------------------------------------
def neurofem_solve(A, b, neurons_per_mesh_point=16, sim_time=2.0):
    """
    Solves Ax = b using the NeuroFEM algorithm exactly as in the notebook.
    Returns the SNN solution vector.
    """
    n = A.shape[0]

    # --------------------------------------------------------------
    # 5.1 Jacobi diagonal preconditioning (not symmetric)
    # --------------------------------------------------------------
    print("  [Preprocess] Jacobi diagonal scaling...")
    diag = A.diagonal()
    if np.any(diag <= 0):
        print(f"  WARNING: {np.sum(diag <= 0)} non‑positive diagonals")
        diag = np.abs(diag)
        diag[diag == 0] = 1.0
    Pinv = sp.sparse.diags(1.0 / diag)
    A_sys = Pinv @ A                     # scaled system matrix
    b_sys = Pinv @ b                     # scaled right‑hand side

    # --------------------------------------------------------------
    # 5.2 Parameters taken directly from the notebook
    # --------------------------------------------------------------
    gamma_pow2 = -6                       # γ = 2^{-6}
    s_gamma = 7 - gamma_pow2
    gamma_norm = (2**gamma_pow2 - 2**-s_gamma)   # ≈ 0.0155029

    omega_s_pow2 = 1
    omega_max = 2**omega_s_pow2 - 2**-(7 - omega_s_pow2)  # ≈ 1.984375

    lambda_d = 8.0        # Hz
    lambda_v = 16.0       # Hz
    ki = 16.0             # integral gain
    kp = 4.0              # proportional gain
    dt = 2**(-12)         # ≈ 0.244 ms
    sigma_v = 0.00225     # noise amplitude
    mu = 0.0              # L2 spike penalty (set to 0)
    nu = 0.0              # L1 spike penalty (set to 0)

    n_timesteps = int(sim_time / dt)

    print(f"  [Parameters] γ = {gamma_norm:.6f}")
    print(f"  [Parameters] ω_max = {omega_max:.6f}")
    print(f"  [Parameters] λ_d = {lambda_d}, λ_v = {lambda_v}")
    print(f"  [Parameters] kp = {kp}, ki = {ki}")
    print(f"  [Parameters] dt = {dt:.6f}, steps = {n_timesteps}")

    # --------------------------------------------------------------
    # 5.3 Compute τ_A (critical scaling factor)
    # --------------------------------------------------------------
    A_max = np.amax(np.abs(A_sys))
    tau_A = (gamma_norm**2) * A_max / omega_max
    print(f"  [Scaling] A_max = {A_max:.6f}, τ_A = {tau_A:.6e}")

    # --------------------------------------------------------------
    # 5.4 Create Γ matrix and weight matrices
    # --------------------------------------------------------------
    gamma = generate_gamma_sparse(n, neurons_per_mesh_point, gamma_norm)
    n_neurons = gamma.shape[1]
    print(f"  [Network] {n_neurons} neurons ({neurons_per_mesh_point} per variable)")

    # NOTE: the notebook passes -A_sys to the weight function
    omega_f, GTAG, threshs, v_reset = create_spiking_fem_network_sparse(
        -A_sys, gamma, lambda_d, lambda_v, mu, nu, tau_A, gamma_norm
    )

    # --------------------------------------------------------------
    # 5.5 Bias = Γ^T * (b / τ_A)
    # --------------------------------------------------------------
    bias = gamma.T @ (b_sys / tau_A)
    print(f"  [Bias] range = [{bias.min():.3f}, {bias.max():.3f}]")
    print(f"  [Threshold] = {threshs[0]:.6f}")

    # --------------------------------------------------------------
    # 5.6 Simulate the spiking network
    # --------------------------------------------------------------
    output, spikes = run_neurofem_simulation(
        GTAG, omega_f, threshs, v_reset, bias,
        gamma, dt, n_timesteps,
        lambda_d, lambda_v, kp, ki, sigma_v
    )
    print(f"  [Simulation] total spikes = {np.sum(spikes)}")

    # --------------------------------------------------------------
    # 5.7 Decode the solution (average over steady state)
    # --------------------------------------------------------------
    burn_in = int(0.2 * n_timesteps)
    x_snn = np.mean(output[:, burn_in:], axis=1)   # solution of the scaled system
    # Jacobi scaling does not change the solution, so x_snn is already the solution of the original system
    return x_snn

# ------------------------------------------------------------------
# 6. MAIN EXECUTION
# ------------------------------------------------------------------
def main():
    model_name = sys.argv[1] if len(sys.argv) > 1 else "Sphere_00"
    print(f"\n{'='*60}")
    print(f"  NeuroFEM (notebook‑exact implementation)")
    print(f"  Problem: {model_name}")
    print(f"{'='*60}")

    # 1. Load data
    print("\n[1] Loading linear system...")
    A, b = load_csr_matrix(model_name)
    n = A.shape[0]
    print(f"   size = {n}×{n}, nnz = {A.nnz}")

    # 2. Solve with NeuroFEM
    print("\n[2] Running NeuroFEM solver...")
    x_snn = neurofem_solve(A, b, neurons_per_mesh_point=16, sim_time=2.0)

    # 3. Compare with exact solution
    print("\n[3] Evaluating accuracy...")
    x_exact = sp.sparse.linalg.spsolve(A, b)

    error = x_exact - x_snn
    mae = np.mean(np.abs(error))
    rmse = np.sqrt(np.mean(error**2))
    rel_err = np.linalg.norm(error) / np.linalg.norm(x_exact) * 100
    corr = np.corrcoef(x_exact, x_snn)[0, 1]

    print(f"   Exact solution: mean = {np.mean(x_exact):.4f}, std = {np.std(x_exact):.4f}")
    print(f"   SNN solution:   mean = {np.mean(x_snn):.4f}, std = {np.std(x_snn):.4f}")
    print(f"   MAE            = {mae:.4f}")
    print(f"   RMSE           = {rmse:.4f}")
    print(f"   Relative error = {rel_err:.2f}%")
    print(f"   Correlation    = {corr:.4f}")

    # 4. Save results
    os.makedirs("data/sol", exist_ok=True)
    out_file_snn = f"data/sol/{model_name}_neurofem_notebook_exact.txt"
    np.savetxt(out_file_snn, x_snn, fmt="%.8f")
    print(f"\n   Solution saved to {out_file_snn}")

    out_file_exact = f"data/sol/{model_name}_exact_solution.txt"
    np.savetxt(out_file_exact, x_exact, fmt="%.8f")
    print(f"\n   Solution saved to {out_file_exact}")

    # 5. Quick plot
    plt.figure(figsize=(12, 4))
    plt.subplot(131)
    plt.plot(x_exact, 'b-', label='Exact', alpha=0.7)
    plt.plot(x_snn, 'r--', label='SNN', alpha=0.7)
    plt.xlabel('Node index')
    plt.ylabel('Solution value')
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.title(f'Solution comparison (corr={corr:.3f})')

    plt.subplot(132)
    plt.scatter(x_exact, x_snn, alpha=0.5, s=20)
    lims = [min(x_exact.min(), x_snn.min()), max(x_exact.max(), x_snn.max())]
    plt.plot(lims, lims, 'k--')
    plt.xlabel('Exact')
    plt.ylabel('SNN')
    plt.axis('equal')
    plt.grid(True, alpha=0.3)
    plt.title('Correlation plot')

    plt.subplot(133)
    plt.hist(error, bins=30, edgecolor='black', alpha=0.7)
    plt.axvline(0, color='r', linestyle='--')
    plt.xlabel('Error')
    plt.ylabel('Frequency')
    plt.grid(True, alpha=0.3)
    plt.title(f'Error distribution (MAE={mae:.3f})')

    plt.tight_layout()
    plt.savefig(f"{model_name}_neurofem_notebook_exact.png", dpi=150)
    print(f"   Plot saved to {model_name}_neurofem_notebook_exact.png")

    print(f"\n{'='*60}")
    print("  Done!")
    print(f"{'='*60}")

if __name__ == "__main__":
    main()