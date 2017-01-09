'''

# A python module for calculating the probability density functions values for 
# location uncertainty assessment
# 
# 
# (C) Zhengguang Zhao, 2016

'''
# !/Users/Uqer/anaconda/bin/python


from numpy import array, square, exp, zeros, floor, subtract, append, sum


def pdftp(N, ngeo, tp, tmp, to, pSigma):
    """
    pdftp() calculates P-wave probability density function values
    N = normalization constant that ensures the integral of the above expression
        over all possible location is equal to one
    ngeo = number of geophones
    tp = computed values of traveltimes by raytracing method 
    tmp = vector of measured arrival times
    to = origin time 
    pSigma = the standard deviations of the measured arrival times

    """

    # what happens here?
    # temp = subtract(tp, tmp) - to
    # temp = square(temp) / (2 * pSigma * pSigma)
    # pdf = N * exp(-sum(temp))

    # realization by loop
    temp = []
    for i in range(ngeo):
        temp.append(tp[i, 0] - tmp[i, 0] - to)
    temp = array(temp)
    temp.shape = (ngeo, 1)
    temp = square(temp) / (2 * pSigma * pSigma)
    n = -sum(temp)
    pdf = N * exp(n)

    return pdf, n
