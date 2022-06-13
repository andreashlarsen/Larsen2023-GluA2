import numpy as np
import matplotlib.pyplot as plt
import os
import time

time_start = time.time()

pepsi='/sansom/s157/bioc1642/Desktop/Scripts/Pepsi-SANS'
data='/sansom/s157/bioc1642/Desktop/prj_GluA2_pH/22/data/Ifilt.dat'
calcfile ='calc.dat'
dir_pep='pepsi'
folder=34 # working directory
FIX_DRHO = 1
PLOT = 0

first     = 1 # start at 1 (not 0)
#last      = 3 # testing/debugging 
last      = 2781 # one more than the actual last because of 0-indexing 
sum_d_rho = 0.0
sum_r0    = 0.0
sum_cst   = 0.0
sum_I0    = 0.0
sum_w     = 0.0

#make directory
os.system('mkdir -p %s' % dir_pep)

# generate calc file (and overwrite old file)
with open(calcfile,'w')as f:
    pass

# import weights
w = np.genfromtxt('weight',usecols=[0],unpack=True)

# fit freely
for i in range(first,last):
    pdb='/sansom/s157/bioc1642/Desktop/prj_GluA2_pH/%d/cg2at/CG2AT%d/PROTEIN/PROTEIN_de_novo_merged.pdb' % (folder,i)
    os.system('%s %s %s -cst -fast -o %s/pepsi_%d.out' % (pepsi,pdb,data,dir_pep,i))

# get average fit params
for i in range(first,last):
    with open('%s/pepsi_%d.log' % (dir_pep,i)) as f:
        
        line = f.readline()
        while line:
            if 'Best d_rho found' in line: 
                line1 = line.split(':')[1]
                d_rho = float(line1.split('e/A')[0])
                print('d_rho = %f' % d_rho)
            if 'Best r0 found' in line: 
                line1 = line.split(':')[1]
                r0    = float(line1.split('A')[0])
                print('r0 = %f' % r0)
            if 'Constant value' in line:
                cst   = float(line.split(':')[1])
            if 'I(0)' in line:
                I0    = float(line.split(':')[1])
            line = f.readline() 
    sum_d_rho += w[i-1]*d_rho
    sum_r0    += w[i-1]*r0
    sum_cst   += w[i-1]*cst
    sum_I0    += w[i-1]*I0
    sum_w     += w[i-1]

mean_d_rho   = sum_d_rho/sum_w
mean_r0      = sum_r0/sum_w
mean_cst     = sum_cst/sum_w
mean_I0      = sum_I0/sum_w

# find default dro and r0
os.system('%s %s %s --dro 1.0 --r0_min_factor 1.0 --r0_max_factor 1.0 --r0_N 1 -fast -o %s/pepsi_default.out' % (pepsi,pdb,data,dir_pep))
with open('%s/pepsi_default.log' % dir_pep, 'r') as f:
    lines = f.readlines()
for line in lines:
    if "Best d_rho found" in line:
        line1=line.split(':')[1]
        default_d_rho=float(line1.split('e/A')[0])
    if "Best r0 found" in line:
        line1=line.split(':')[1]
        default_r0=float(line1.split('A')[0])

print('\n\n\n\n\n\n\n')
print('************\n')
print('SECOND round\n')
print('mean_d_rho = %f' % mean_d_rho)
print('mean_r0 = %f' % mean_r0)
print('mean_cst = %f' % mean_cst)
print('mean_I0 = %f' % mean_I0)
print('default_d_rho = %f' % default_d_rho)
print('default_r0 = %f' % default_r0)
print('************\n')
print('\n\n\n\n\n\n\n')

for i in range(first,last):
    pdb='/sansom/s157/bioc1642/Desktop/prj_GluA2_pH/%d/cg2at/CG2AT%d/PROTEIN/PROTEIN_de_novo_merged.pdb' % (folder,i)
    if FIX_DRHO:
        os.system('%s %s %s --dro %f --r0_min_factor %f --r0_max_factor %f --r0_N 1 -cst --cstFactor %f --I0 %f -fast -o %s/pepsi_const_%d.out' % (pepsi,pdb,data,5.0,mean_r0/default_r0,mean_r0/default_r0,mean_cst,mean_I0,dir_pep,i))
    else:
        os.system('%s %s %s --dro %f --r0_min_factor %f --r0_max_factor %f --r0_N 1 -cst --cstFactor %f --I0 %f -fast -o %s/pepsi_const_%d.out' % (pepsi,pdb,data,mean_d_rho/default_d_rho,mean_r0/default_r0,mean_r0/default_r0,mean_cst,mean_I0,dir_pep,i))

# write to calcfile and calc mean_Ifit
sum_Ifit = 0.0
for i in range(first,last):
    Ifit = np.genfromtxt('%s/pepsi_const_%d.out' % (dir_pep,i),skip_header=6,usecols=[3],unpack=True)
    with open(calcfile,'a') as f:
        f.write('frame%d ' %i)
        for j in range(len(Ifit)):
            f.write('%f ' % Ifit[j])
        f.write('\n')
    sum_Ifit += w[i-1]*Ifit
mean_Ifit = Ifit/sum_w

q,Idat,sigma = np.genfromtxt('%s/pepsi_const_%d.out' % (dir_pep,i),skip_header=6,usecols=[0,1,2],unpack=True)

R = (Ifit-Idat)/sigma
chi2r = np.sum(R**2)/len(Idat)

# plot mean fit before BME
if PLOT:
    plt.errorbar(q,Idat,yerr=sigma,linestyle='none',marker='.',color='red',label='SANS, GluA2 AMPA pH 5.5')
    plt.plot(q,Ifit,color='black',label='Fit, before BME, chi2r=%1.2f' % chi2r)
    plt.xscale('log')
    plt.yscale('log')
    plt.legend()
    plt.show()

## timing
time_elapsed = time.time()-time_start
print('time for pepsi.py: %1.2f s' % time_elapsed)

