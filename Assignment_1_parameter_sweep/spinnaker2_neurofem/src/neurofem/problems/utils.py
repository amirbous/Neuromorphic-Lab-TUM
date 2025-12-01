import numpy as np
import scipy.sparse as sp

def float_to_signed_sparse(matrix: sp.spmatrix | np.ndarray, x_bits: int = 21, scale: float | None = None) -> tuple[sp.spmatrix, float]:
    """
    Quantizize a sparse matrix to signed integers with specific bit-width.

    The input floating point matrix is converted to a sparse integer matrix using uniform quantization: quantized_value = round(original_value / scale)

    The range of representable integers is determined by the number of bits specified (x_bits). One bit is reserved for the sign.
    For example, for x_bits=6, the representable range is from -32 to 31.

    :param matrix: The input sparse matrix.
    :param x_bits: Number of bits for the quantized integers (including sign bit).
    :param scale: Optional scale factor for quantization. The matrix is divied by this value before rounding. Will be computed if not provided.
    :return: The quantized CSR matrix of type int32 and the scale factor.
    """
    if not sp.issparse(matrix):
        matrix = sp.csr_matrix(matrix)

    x_bits -= 1 # one bit for sign

    matrix = matrix.astype(np.float32).tocoo()
    max_int = 2**x_bits - 1
    min_int = -2**x_bits

    # compute scale if not given
    if scale is None:
        max_val = np.max(np.abs(matrix.data))
        scale = max_val / max_int if max_val != 0 else 1.0

    # efficient quantitization in sparse form
    int_data = np.empty_like(matrix.data, dtype=np.int32)
    for i in range(len(matrix.data)):
        val = matrix.data[i] / scale
        int_data[i] = int(np.clip(np.round(val), min_int, max_int))

    int_matrix = sp.coo_matrix((int_data, (matrix.row, matrix.col)), shape=matrix.shape, dtype=np.int32).tocsr()
    return int_matrix, scale
