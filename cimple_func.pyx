def hamming (str x1, str x2, int m):
    cdef j = 0
    for i in range(len(x1)):
        if x1[i] != x2[i]:
            j = j + 1
            if j > m:
                return (False)
    return (True)