import sys,os
import numpy as np
sys.path.append('/sansom/s157/bioc1642/Desktop/prj_GluA2_pH/22/BME/BME2-master2/')
import BME
import matplotlib.pyplot as plt

PLOT = 0

# define input file names
exp_file = '/sansom/s157/bioc1642/Desktop/prj_GluA2_pH/22/data/Ifilt.dat'
calc_file = '../calc.dat'
weight_file = '../weight'
initial_weights = np.loadtxt(weight_file)

# BME settings
prefix = 'gluA2'
maxite = 50
ftol   = 0.001

# initialize. A name must be specified
print('INITIALIZING BME')
rew = BME.Reweight(prefix,w0=initial_weights)

# load the experimental and calculated datasets
rew.load(exp_file,calc_file,fit="scale+offset")

#rew.ibme(theta=1000,iterations=50,ftol=0.00001)

#results = rew.fit(theta=100)

#print("CHI2  original: %6.2f" % results[0])
#print("CHI2 optimized: %6.2f" % results[1])

n_theta=200
thetas_exp = np.linspace(-1,4,num=n_theta)
thetas = 10**thetas_exp
#fig,ax = plt.subplots(figsize=(20, 10))

print('BME THETA loop')
print('thetas:')
print(thetas)
f = open('chi_phi.dat','w')
f.write('theta   chi2r_after, chi2r_before, phi\n')
chi2_ratio = []
phis = []
chi2 = []
count = 0
for t in thetas:
    count += 1
    print('theta = %f' % t)
    print('iteration %d/%d' % (count,n_theta))
    chi2_before,chi2_after,phi = rew.ibme(theta=t,iterations=maxite,ftol=ftol)
    phis.append(phi)
    chi2.append(chi2_after)
    chi2_ratio.append(chi2_after/chi2_before)
    f.write('%f %e %e %e\n' % (t,chi2_after,chi2_before,phi))
f.close()

if PLOT:
    plt.plot(thetas,phis,"-o",label="Phi",c="k")
    plt.plot(thetas,chi2_ratio,"-o",label="Chi2 reduction",c="r")
    plt.plot(thetas,chi2,"-o",label="Chi2",c="green")
    plt.plot([0.01,10000],[1.0,1.0],"--",c="grey")
    plt.plot([0.01,10000],[1.6,1.6],"--",c="black",label='Chi2 from BIFT')
    plt.legend()
    plt.xscale('log')
    plt.xlabel("Theta")
    plt.legend()
    plt.show()
    plt.savefig('BME_gluA2.pdf',format='pdf')

os.system('rm %s_*' % prefix)

for theta_opt in [1]:
#for theta_opt in [5,10,15,20,25,30,35,40,45,50,60,70,80,90,100,150,200,300,400,500,600,700,800,900,1000,10000]:
#for theta_opt in np.arange(0.1,10.0,0.2):   
    print('theta_opt = %f' % theta_opt)
    chi2_before,chi2_after,phi = rew.ibme(theta=theta_opt,iterations=maxite,ftol=ftol)
    print('phi = %e\n\n' % phi)
    os.system('rm %s*log' % prefix)
    os.system('mv %s*.calc.dat theta%1.1f.calc.dat' % (prefix,theta_opt))
    os.system('mv %s*.weights.dat theta%1.1f.weights.dat' % (prefix,theta_opt))

