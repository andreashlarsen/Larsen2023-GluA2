import numpy as np
import matplotlib.pyplot as plt

for i in range(270,279):
    
    ## import data
    filename = 'pepsi_%d.out' %i
    filepath = 'pepsi/%s' % filename
    q,I,dI,Ifit = np.genfromtxt(filepath,skip_header=6,usecols=[0,1,2,3],unpack=True)

    ## plot data and fit
    plt.errorbar(q,I,yerr=dI,linestyle='none',marker='.',color='red',label=filepath)
    plt.plot(q,Ifit,color='black',label='fit')

    ## figure settings
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('q')
    plt.ylabel('I')
    plt.legend(frameon=False)
    plt.show()

