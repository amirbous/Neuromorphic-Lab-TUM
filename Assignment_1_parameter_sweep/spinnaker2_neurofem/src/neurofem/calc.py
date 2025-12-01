import numpy as np

def calc_mse(matrix: np.ndarray, rhs: np.ndarray, solution: np.ndarray) -> float:
    """
    Calculate the mean squared error (MSE) between Ax and b.

    :param matrix: The system matrix A.
    :param rhs: The right-hand side vector b.
    :param solution: The estimated solution vector x.
    :return: The MSE value.
    """
    residual = matrix @ solution - rhs
    mse = np.mean(residual ** 2)
    return mse

def calc_rmse(matrix: np.ndarray, rhs: np.ndarray, solution: np.ndarray) -> float:
    """
    Calculate the root mean squared error (RMSE) between Ax and b.

    :param matrix: The system matrix A.
    :param rhs: The right-hand side vector b.
    :param solution: The estimated solution vector x.
    :return: The RMSE value.
    """
    mse = calc_mse(matrix, rhs, solution)
    rmse = np.sqrt(mse)
    return rmse

def calc_relative_residual(matrix: np.ndarray, rhs: np.ndarray, solution: np.ndarray) -> float:
    """
    Calculate the relative residual ||Ax - b|| / ||b||.

    :param matrix: The system matrix A.
    :param rhs: The right-hand side vector b.
    :param solution: The estimated solution vector x.
    :return: The relative residual value.
    """
    residual = matrix @ solution - rhs
    relative_residual = np.linalg.norm(residual) / np.linalg.norm(rhs)
    return relative_residual
