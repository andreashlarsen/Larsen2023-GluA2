import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

PLOT = 0

## first and last data points used in Guinier approximation
first = 2
last = 26
Rg0  = 55

f = open('rg_pepsi.dat','w')

sum_Rg = 0.0
count = 0
min_Rg = 1e8
max_Rg = 0.0
for i in range(1,2781):
#for i in range(275,279):
    q,Ipepsi = np.genfromtxt('pepsi/pepsi_%d.out' %i,skip_header=6,usecols=[0,3],unpack=True)

    ## get Rg by Guinier approximation
    lnI = np.log(Ipepsi)
    q2  = q**2
    maxq = q[last]
    def func(q2,lnI0,Rg2,p0=[lnI[0],Rg0]):
        return lnI0-q2*Rg2/3.0
    popt,popv = curve_fit(func,q2[first:last],lnI[first:last])
    Rg = np.sqrt(popt[1])
    maxqRg = maxq*Rg
    print('i= %d, Rg = %1.1f, max(q)*Rg = %1.1f' % (i,Rg,maxqRg))
    f.write('%f\n' % Rg)

    lnIfit=func(q2,*popt)
    Ifit = np.exp(lnIfit)

    sum_Rg += Rg
    count += 1
    if Rg < min_Rg:
        min_Rg = Rg
    if Rg > max_Rg:
        max_Rg = Rg
    
    ## plot the Guinier fit
    if PLOT:
        fig,(ax0,ax1) = plt.subplots(2,1,gridspec_kw = {'height_ratios':[3,1]})
        ax0.plot(q2,lnI,linestyle='none',marker='.',label='Pepsi fit')
        ax0.plot(q2[first:last],lnI[first:last],linestyle='none',marker='o',label='Included in Guinier fit')
        ax0.plot(q2,lnIfit,color='darkgreen',label='Guinier fit')
        ax0.set_xlabel(r'$q^2$')
        ax0.set_ylabel('ln(I)')

        ax0.set_xlim(0,0.002)
        ax0.set_ylim(-6,-3)
        ax0.legend()

        R = lnI-lnIfit
        maxR = np.ceil(np.amax(abs(R*10)[first:last]))/10
        print(maxR)

        ax1.plot(q2,R,linestyle='none',marker='.')
        ax1.plot(q2[first:last],R[first:last],linestyle='none',marker='o')
        ax1.plot(q2,q2*0,color='darkgreen')
        ax1.set_xlim(0,0.002)
        ax1.set_ylim(-maxR,maxR)

        plt.show()
    
f.close()

mean_Rg = sum_Rg/count
print('mean Rg = %1.1f' % mean_Rg)
print('min Rg = %1.1f' % min_Rg)
print('max Rg = %1.1f' % max_Rg)

