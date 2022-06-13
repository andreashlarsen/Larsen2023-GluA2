import sys,os
import numpy as np
import matplotlib.pyplot as plt
from rebin import rebin

FITS = 0
RG = 0
RGLOG = 0

#Rg_4u2p=54.4
Rg_exp=60.57
theta_opt=1.0
n_frames=2780 # number of frames

print('Theta_optimal = %1.1f' % theta_opt)

w0 = np.genfromtxt('../weight',usecols=[0],unpack=True)
sum_w0 = np.sum(w0)
w0 /= sum_w0
frame,w = np.genfromtxt('theta%1.1f.weights.dat' % theta_opt,usecols=[0,1],unpack=True)

# find non-zero weights
idx = np.where(w>0.0)

## import pepsi Rg (by Guinier fit)
rg_pepsi = np.genfromtxt('../rg_pepsi.dat',usecols=[0],unpack=True)
print('min Rg pepsi = %1.2f' % np.amin(rg_pepsi))
print('max Rg pepsi = %1.2f' % np.amax(rg_pepsi))

## import dist
#time,dist = np.genfromtxt('../DIST',comments='#',usecols=[0,1],unpack=True)
dist = np.genfromtxt('../distance',usecols=[0],unpack=True)

## sort out Rg outliers
threshold = 1000
idx_rg = np.argwhere(rg_pepsi<threshold)
print('max Rg pepsi (threshold) = %1.2f' % np.amax(rg_pepsi[idx_rg]))

# Rg with different weights
Rg_w0 = np.sum(rg_pepsi[idx_rg]*w0[idx_rg])
Rg_w = np.sum(rg_pepsi[idx_rg]*w[idx_rg])
print('Rg*w0 = %1.1f, Rg*w = %1.1f, Rg_exp = %1.1f+/-0.7' % (Rg_w0,Rg_w,Rg_exp))

row=3
col=1
p1  = plt.subplot(row,col,1)
p2a = plt.subplot(row,col,2)
p2b = p2a.twinx()
p3  = plt.subplot(row,col,3)

## Weight vs Frame
bar_width = 20
opacity = 0.5
p1.bar(frame,w0,width=bar_width,alpha=opacity,color='red',label='w0')
p1.bar(frame+bar_width,w,width=bar_width,alpha=opacity,color='blue',label='w') 
p1.set_ylabel('Weight')
p1.set_xlabel('Frame')
p1.set_xlim(0,n_frames)
p1.legend(loc='upper right',frameon=False)
#p1.set_yscale('log')

# Rg vs Frame
p2a.plot(frame[idx_rg],rg_pepsi[idx_rg],color='darkcyan')
lw=2
p2a.axhspan(Rg_exp-0.7,Rg_exp+0.7,alpha=0.5,color='grey',label=r'$R_{g,exp}$')
p2a.plot([1,n_frames],[Rg_w0,Rg_w0],linestyle='--',linewidth=lw,color='red',label=r'$R_{g,w0}$')
p2a.plot([1,n_frames],[Rg_w,Rg_w],linestyle='--',linewidth=lw,color='blue',label=r'$R_{g,w}$')
p2a.set_ylabel(r'$R_g$ [$\AA$]',color='darkcyan')
p2a.set_xlabel('Frame')
p2a.set_xlim(0,n_frames)
p2a.legend(loc='upper right',frameon=False)

# NTD vs frame (same axes)
p2b.plot(frame,dist,color='black',label='NTD dist')
p2b.set_ylabel(r'NTD dist [$\AA$]',color='black')
p2b.set_ylim(0,25)

# Frequency vs Rg
bins = 50
alpha = 0.5
y = [1e-80,1e20]
p3.hist(rg_pepsi[idx_rg],weights=w0[idx_rg],bins=bins,color='red',alpha=alpha)#,label='w0')
p3.hist(rg_pepsi[idx_rg],weights=w[idx_rg],bins=bins,color='blue',alpha=alpha)#,label='w')
p3.plot([Rg_w0,Rg_w0],y,linestyle='--',linewidth=lw,color='red',label=r'$R_{g,w0}$')
p3.plot([Rg_w,Rg_w],y,linestyle='--',linewidth=lw,color='blue',label=r'$R_{g,w}$')
plt.axvspan(Rg_exp-0.7,Rg_exp+0.7,alpha=0.5,color='grey',label=r'$R_{g,exp}$')
p3.legend(loc='upper right',frameon=False)
p3.set_ylabel('Frequency')
p3.set_xlabel(r'$R_g$ [$\AA$]')
p3.set_yscale('log')
p3.set_ylim(y)

#plt.tight_layout()
plt.savefig('BME_weight.pdf',format='pdf')
plt.show()

if FITS:
    # define input file names
    exp_file = '/sansom/s157/bioc1642/Desktop/prj_GluA2_pH/22/data/Ifilt.dat'
    q,I,sigma = np.genfromtxt(exp_file,skip_header=1,usecols=[0,1,2],unpack=True)
    M = len(q)

    #calculate fits before and after BME
    Iw = np.zeros(M)
    for i in range(M):
        Ifit = np.genfromtxt('theta%1.1f.calc.dat' % theta_opt,usecols=[i+1],unpack=True)
        Iw[i] = np.sum(w[idx_rg]*Ifit[idx_rg])

    Iw0 = np.zeros(M)
    for i in range(M):
        Ifit = np.genfromtxt('theta%1.1f.calc.dat' % 1e4,usecols=[i+1],unpack=True)
        Iw0[i] = np.sum(w0[idx_rg]*Ifit[idx_rg])
    Rw = (Iw-I)/sigma
    Rw0 = (Iw0-I)/sigma
    chi2r_w = np.sum(Rw**2)/(M-2)
    chi2r_w0 = np.sum(Rw0**2)/(M-2)

    f,(p1,p2) = plt.subplots(2,1,gridspec_kw={'height_ratios':[4,1]})

    p1.errorbar(q,I,yerr=sigma,linestyle='none',marker='.',color='grey',zorder=0)
    p1.plot(q,Iw0,color='red',zorder=1,label=r'Fit before BME reweighting, $\chi^2_r = %1.1f$' % chi2r_w0)
    p1.plot(q,Iw,color='blue',zorder=2,label=r'Fit after BME reweighting, $\chi^2_r = %1.1f$' % chi2r_w)

    p1.set_xscale('log')
    p1.set_yscale('log')
    p1.set_ylabel(r'$I(q)$ [cm$^{-1}$]')
    p1.set_xticks([])
    p1.legend(frameon=False)

    p2.plot(q,np.zeros(M),linestyle='none',marker='.',color='grey')
    p2.plot(q,Rw0,color='red')
    p2.plot(q,Rw,color='blue')
    p2.set_xscale('log')
    p2.set_xlabel(r'$q$ [$\AA^{-1}$]')
    p2.set_ylabel(r'$\Delta I/\sigma$')
    Rmax = np.ceil(np.amax(Rw0))
    p2.set_ylim(-Rmax,Rmax)
    p2.set_yticks([-Rmax,0,Rmax])

    plt.show()

if RG:
    bins = 50
    alpha = 0.5
    plt.hist(rg_pepsi[idx_rg],weights=w0[idx_rg],bins=bins,color='red',alpha=alpha)#,label='w0')
    plt.hist(rg_pepsi[idx_rg],weights=w[idx_rg],bins=bins,color='blue',alpha=alpha)#,label='w')    
    plt.axvspan(Rg_exp-0.7,Rg_exp+0.7,alpha=0.5,color='grey',label=r'$R_{g,exp}$')
    #plt.plot([Rg_exp,Rg_exp],[0,ymax],linestyle='--',linewidth=lw,color='black',label=r'$R_{g,exp}$')
    plt.plot([Rg_w0,Rg_w0],[0,ymax],linestyle='--',linewidth=lw,color='red',label=r'$R_{g,w0}$')
    plt.plot([Rg_w,Rg_w],[0,ymax],linestyle='--',linewidth=lw,color='blue',label=r'$R_{g,w}$')
    plt.legend(loc='upper right',frameon=False)
    plt.ylabel('Frequency')
    plt.xlabel(r'$R_g$ [$\AA$]')
    plt.ylim(0,ymax)
    plt.xlim(50,65)

    plt.show()

if RGLOG:
    bins = 50
    alpha = 0.5
    ymin = 1e-80
    ymax = 1e10
    plt.hist(rg_pepsi[idx_rg],weights=w0[idx_rg],bins=bins,color='red',alpha=alpha)#,label='w0')
    plt.hist(rg_pepsi[idx_rg],weights=w[idx_rg],bins=bins,color='blue',alpha=alpha)#,label='w')    
    plt.axvspan(Rg_exp-0.7,Rg_exp+0.7,alpha=0.5,color='grey',label=r'$R_{g,exp}$')
    plt.plot([Rg_w0,Rg_w0],[0,ymax],linestyle='--',linewidth=lw,color='red',label=r'$R_{g,w0}$')
    plt.plot([Rg_w,Rg_w],[0,ymax],linestyle='--',linewidth=lw,color='blue',label=r'$R_{g,w}$')
    plt.legend(loc='upper right',frameon=False)
    plt.ylabel('Frequency')
    plt.xlabel(r'$R_g$ [$\AA$]')
    plt.yscale('log')
    plt.ylim(ymin,ymax)
    plt.xlim(50,65)

    plt.show()
