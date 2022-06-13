import sys,os
import numpy as np
import matplotlib.pyplot as plt

theta,chi2_after,chi2_before,phi = np.genfromtxt('chi_phi.dat',skip_header=1,usecols=[0,1,2,3],unpack=True)
S = np.log(phi)

theta_opt=1.0
min_S = -20
max_S = 120
min_theta = 0.1
max_theta = 1000
min_chi2 = 0
max_chi2 = 9

idx = np.where(theta>=theta_opt)[0][0]
#print(idx)
print('theta_opt,-S_opt,Chi2r_opt = %1.1f,%1.2f,%1.1f' % (theta[idx],-S[idx],chi2_after[idx]))
p1 = plt.subplot(1,2,1)
p2 = p1.twinx()
p3 = plt.subplot(1,2,2)

p1.plot(theta,chi2_after,label=r'$\chi^2_r$',color='grey')
p1.plot([min_theta,max_theta],[1.6,1.6],'--',color='grey',label='Chi2 from BIFT')
p1.set_xscale('log')
p1.set_xlabel(r'$\theta$')
p1.set_ylabel(r'$\chi^2_r$',color='grey')
p1.set_ylim(min_chi2,max_chi2)
p1.set_xlim(min_theta,max_theta)
p1.tick_params(axis='y',colors='grey')
p1.plot([theta_opt,theta_opt],[min_chi2,max_chi2],'-',color='black')

p2.plot(theta,-S,label=r'$S_\mathrm{REL}$',color='green')
p2.set_ylabel(r'$-S_\mathrm{REL}$',color='green')
p2.tick_params(axis='y', colors='green')
p2.plot([min_theta,max_theta],[0.0,0.0],'--',color='green',label='S_REL = 0')
p2.set_ylim(min_S,max_S)

p3.plot(-S,chi2_after,color='black')
p3.plot([min_S,max_S],[1.6,1.6],'--',color='grey',label='Chi2 from BIFT')
p3.plot([0,0],[min_chi2,max_chi2],'--',color='green',label='S_REL = 0')
p3.plot(-S[idx],chi2_after[idx],linestyle='none',marker='o',color='black')
p3.set_xlabel(r'-$S_\mathrm{REL}$',color='green')
p3.tick_params(axis='x',colors='green')
p3.tick_params(axis='y',colors='grey')
p3.set_ylabel(r'$\chi^2_r$',color='grey')
p3.set_ylim(min_chi2,max_chi2)
p3.set_xlim(min_S,110)

plt.tight_layout()
plt.savefig('BME_gluA2.pdf',format='pdf')
plt.show()

