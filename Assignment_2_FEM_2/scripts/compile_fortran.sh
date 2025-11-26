#!/usr/bin/env bash
set -euo pipefail

# compile_fortran.sh
# Purpose: Compile the Fortran sources in `src/` into an executable named `main`.
# Usage: ./compile_fortran.sh
# Notes:
#  - Uses `gfortran` with optimization level `-O2` and Fortran 2008 standard.
#  - Script determines the project root relative to its own location and runs
#    the compiler in the `src/` directory so file paths are predictable.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
SRC_DIR="$ROOT_DIR/src"

echo "[compile_fortran] Project root: $ROOT_DIR"
echo "[compile_fortran] Source directory: $SRC_DIR"

if ! command -v gfortran >/dev/null 2>&1; then
  echo "Error: gfortran not found in PATH. Install GNU Fortran and retry." >&2
  exit 1
fi

cd "$SRC_DIR"

echo "[compile_fortran] Running gfortran..."
gfortran -O2 -std=f2008 mod_types.f90 mod_rhs.f90 mod_mesh.f90 mod_assemble.f90 mod_solver.f90 mod_io.f90 mod_compute_error.f90 main.f90 -o main

if [ $? -eq 0 ]; then
  echo "[compile_fortran] Compilation succeeded. Executable: $SRC_DIR/main"
else
  echo "[compile_fortran] Compilation failed." >&2
  exit 1
fi

exit 0
