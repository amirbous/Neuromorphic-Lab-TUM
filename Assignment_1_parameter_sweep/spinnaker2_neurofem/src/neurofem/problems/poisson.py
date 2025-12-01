import numpy as np
import skfem as fem
from skfem.models import poisson

def generate_poisson_unitdisk_variable_f(nrefs: int = 5):
    """
    Assemble the linear system Ax=b for a Poisson problem on the unit disk.
    The continuous problem is:
        -Δu = f in unit disk, u=0 on boundary
    with:
        f(x,y) = 12 - 60 * (x - 0.25)^2 - 60 * (y + 0.13)^2

    :param nrefs: Number of mesh refinements.
    :return: A tuple (matrix, rhs, basis) where
        - matrix is the assembled stiffness matrix A (CSR format)
        - rhs is the assembled right-hand side vector b
        - basis is the finite element basis object
    """
    mesh = fem.MeshTri().init_circle(nrefs)
    basis = fem.CellBasis(mesh, fem.ElementTriP1())

    def f_fun(v ,w):
        x, y = w.x
        f_val = -(12 - 60 * (x - 0.25) ** 2 - 60 * (y + 0.13) ** 2)
        return f_val * v

    matrix = poisson.laplace.assemble(basis)
    rhs = fem.asm(fem.LinearForm(f_fun), basis)

    # Assemble dirichlet boundary conditions (u=0 on boundary)
    D = basis.get_dofs()
    matrix, rhs = fem.enforce(matrix, rhs,D=D)

    return matrix.tocsr(), rhs, basis
