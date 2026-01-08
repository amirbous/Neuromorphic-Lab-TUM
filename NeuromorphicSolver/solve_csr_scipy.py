import sys
import os
import numpy as np
import scipy as sp
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import spsolve

def main():


    problem_name = sys.argv[1] if len(sys.argv) > 1 else "Sphere_00"
    
    mtx_filename = f"{problem_name}_mtx.txt"
    rhs_filename = f"{problem_name}_RHS.txt"
    sol_filename = f"{problem_name}_sol.txt"


    if not os.path.exists(mtx_filename):
        print(f"Error: {mtx_filename} not found.")
        sys.exit(1)

    print(f"Loading matrix from {mtx_filename}...")
    try:
    
        raw_data = np.loadtxt(mtx_filename, skiprows=1, dtype=float)
        
        
        with open(mtx_filename, 'r') as f:
            header = f.readline().split()
            n_rows = int(header[0])
            n_cols = int(header[1])


        rows = raw_data[:, 0].astype(int) - 1
        cols = raw_data[:, 1].astype(int) - 1
        values = raw_data[:, 2]


        A = sp.sparse.csr_matrix((values, (rows, cols)), shape=(n_rows, n_cols))


        initial_nnz = A.nnz
        A.eliminate_zeros()
        print(f"Shape: {A.shape}")
        print(f"Non-zeros: {A.nnz}")

    except Exception as e:
        print(f"Failed to load matrix: {e}")
        sys.exit(1)


    try:


        b = np.loadtxt(rhs_filename, skiprows=1)
        
        if A.shape[0] != b.shape[0]:
            print(f"Error: Matrix rows ({A.shape[0]}) != RHS length ({b.shape[0]})")
            sys.exit(1)

    except Exception as e:
        print(f"Failed to load RHS: {e}")
        sys.exit(1)


    try:


        x = spsolve(A, b)
        
        if np.any(np.isnan(x)) or np.any(np.isinf(x)):
            print("Error: Solver returned NaN or Inf. Matrix might be singular.")
            sys.exit(1)

    except Exception as e:
        print(f"Solver failed: {e}")
        sys.exit(1)


    print(f"Writing solution to {sol_filename}...")
    try:
        with open(sol_filename, 'w') as f:
            f.write(f"{len(x)}\n")
            

            for val in x:
                f.write(f"{val:.16e}\n")


    except Exception as e:
        print(f"Failed to write solution: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()
