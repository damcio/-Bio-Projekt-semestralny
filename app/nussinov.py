import numpy as np
from app.helpers.inputParser import base_pair_tuple


def w(el1, el2):
    """
    Returns weight based on residue pair
    """
    pair = el1 + el2
    return int(pair == 'AU' or
               pair == 'UA' or
               pair == 'GC' or
               pair == 'CG')


def nussinov_algorithm(rna_strand):
    """
    Finds strands secondary structure using Nussinov's dynamic programming algorithm
    :param rna_strand: Text representation of RNA strand
    :return: List of found base pairs
    """
    matrix = np.zeros(shape=(len(rna_strand), len(rna_strand)), dtype=int)
    for i in range(len(rna_strand)):
        for j in range(len(rna_strand)):
            if i - 1 > j:
                matrix[i, j] = -1

    for i in range(1, len(rna_strand)):
        for j in range(i, len(rna_strand)):
            v1 = matrix[j - i + 1, j - 1] + w(rna_strand[j - i], rna_strand[j])
            v2 = matrix[j - i + 1, j]
            v3 = matrix[j - i, j - 1]
            v4 = max([matrix[j - i, k] + matrix[k + 1, j] for k in range(j - i + 1, j)] + [-1])
            matrix[j - i, j] = max(v1, v2, v3, v4)
    solution = []

    def traceback(i, j):
        if i >= j:
            return solution
        elif matrix[i, j] == matrix[i + 1][j]:
            traceback(i + 1, j)
        elif matrix[i][j] == matrix[i][j - 1]:
            traceback(i, j - 1)
        elif matrix[i][j] == matrix[i + 1][j - 1] + w(rna_strand[i], rna_strand[j]):
            solution.append((i, j))
            traceback(i + 1, j - 1)
        else:
            val, k = max([(matrix[i][k] + matrix[k + 1][j], k) for k in range(i + 1, j)] + [(-1, -1)])
            assert matrix[i][j] == val
            traceback(i, k)
            traceback(k + 1, j)

    traceback(0, len(rna_strand) - 1)

    return [(base_pair_tuple('A', pair[0]+1), base_pair_tuple('A', pair[1]+1)) for pair in solution]
