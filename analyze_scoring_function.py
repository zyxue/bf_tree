def f(a, c, k=20, l=150):
    d = float(l - k)

    s = 0
    for i in xrange(c):
        for j in xrange(a):
    j = 0
    for kmer in U.kmerize(read, k):
        if kmer in bf:
            j += 1
            s += (1 - 1 / (j + 1)) / d
        else:
            j = 0
    return s
                                        
