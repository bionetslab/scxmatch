from math import comb, factorial, log, exp
from itertools import chain


def _cross_match_count(Z, matching, test_group):
    print("counting cross matches.")
    pairs = [(Z.iloc[i], Z.iloc[j]) for (i, j) in matching]
    filtered_pairs = [pair for pair in pairs if (pair[0] == test_group) ^ (pair[1] == test_group)] # cross-match pairs contain test group exactly once
    a1 = len(filtered_pairs)
    return a1


def _get_p_value(a1, n, N, I):
    p_value = 0
    for A1 in range(a1 + 1):  # For all A1 <= a1
        A2 = (n - A1) / 2 
        A0 = I - (n + A1) / 2 

        if int(A0) != A0:
            continue
        if int(A2) != A2:
            continue 
        if A0 < 0 or A2 < 0:
            continue  

        A0 = int(A0)
        A2 = int(A2)
        
        log_numerator = A1 * log(2) + log(factorial(I))
        log_denominator = log(comb(N, n)) + log(factorial(A0)) + log(factorial(A1)) + log(factorial(A2))
        p_value += exp(log_numerator - log_denominator)

    return p_value

    
def _get_z_score(a1, n, N):
    m = N - n
    E = n * m / (N - 1) # Eq. 3 in Rosenbaum paper
    var = 2 * n * (n - 1) * m * (m - 1) / ((N - 3) * (N - 1)**2)
    z = (a1 - E) / np.sqrt(var)
    return z


def _get_relative_support(N, Z):
    return N / len(Z)
    

def _rosenbaum_test(Z, matching, test_group):
    used_elements = list(chain.from_iterable(matching))
    n = sum(1 for el in used_elements if Z.iloc[el] == test_group)
    N = len(matching) * 2
    I = len(matching)

    a1 = _cross_match_count(Z, matching, test_group)
    
    p_value = _get_p_value(a1, n, N, I)
    z_score = _get_z_score(a1, n, N)
    relative_support = _get_relative_support(N, Z)
    return p_value, z_score, relative_support

