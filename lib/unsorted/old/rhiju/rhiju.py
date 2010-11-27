from math import *

def hello():
    print 'Hello World'


def erfcc(x):
    # This function is a copy of the erfcc from ROSETTA, which
    # is used to calculate maxsub approximate maxsub expectation values!
    z = abs(x)
    t = 1.0/(1.0+0.5*z)
    erfval = t*exp(-z*z-1.26551223+t*(1.00002368+t*(.37409196+
	 t*(.09678418+t*(-.18628806+t*(.27886807+t*(-1.13520398+t*(1.48851587+
	 t*(-.82215223+t*.17087277)))))))))
    if ( x < 0.0 ):
        erfval = 2.0 - erfval
    return erfval



def maxlgE(maxsub, proteinlength):
    #Copied from Rosetta ... estimates expectation
    # values from a Gaussian approximation for MAXSUB distribution.
    # I can't seem to find the paper where these magic values came from, but
    # suspect that Angel Ortiz used this in PCONS.

    psi = (100.0*maxsub)/proteinlength

    am = 747.29 * proteinlength**(-0.7971)
    as = 124.99 * proteinlength**(-0.6882)

    zscore = (psi - am) / as

    evalue = 0.5 * erfcc( zscore / sqrt(2))

    znew = 0.730* ((1.2825755*zscore) + 0.5772);
    evalue = 1.0 - exp( - exp(-znew))
    if (evalue < 2.650E-14) : evalue=2.650E-14

    logeval = log(evalue)

    return logeval

