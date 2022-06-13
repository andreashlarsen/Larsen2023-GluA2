import numpy as np

def rebin(x,Nbin):
    
    Nx = len(x)
    j = 0
    Nx_RB = int(np.ceil(Nx/Nbin))
    x_RB = np.zeros(Nx_RB)
    for j in range(Nx_RB):
        rest = Nx-j*Nbin
        if rest > Nbin:
            imax = Nbin
        else:
            imax = rest
        for i in range(imax):
            x_RB[j] += x[j+i]
        j += 1

    return x_RB
