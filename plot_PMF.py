import numpy as np
import matplotlib.pyplot as plt

PLOT_ALL = 0

elements =[\
        #['7a/axon_umbrella_prot7a_glua2_ph','7: NE2, default (neutral)','grey',7],\
        #['7b/axon_umbrella_prot7b_glua2_ph','7b: NE2, default (neutral)','grey',7],\
        #['7c/axon_umbrella_prot7c_glua2_ph','7b: NE2, default (neutral)','grey',7],\
        #['7d/axon_umbrella_prot7d_glua2_ph','7d: NE2, default (neutral)','grey',7],\
        #['7e/axon_umbrella_prot7e_glua2_ph','7e: NE2, default (neutral)','grey',7],\
        #['11/axon_umbrella_prot11_glua2_ph','11: NE2, 20 ns','black',11],\
        #['6a/axon_umbrella_prot6a_glua2_ph','6a: D1, (alternative neutral)','orange',6],\
        #['6b/axon_umbrella_prot6b_glua2_ph','6b: D1, (alternative neutral)','orange',6],\
        #['6c/axon_umbrella_prot6c_glua2_ph','6c: D1, (alternative neutral)','orange',6],\
        #['8/axon_umbrella_prot8_glua2_ph','8: NE2+D1 (acidic)','green',8],\
        #['8a/axon_umbrella_prot8a_glua2_ph','8a: NE2+D1 (acidic)','green',8],\
        #['8b/axon_umbrella_prot8b_glua2_ph','8b: NE2+D1 (acidic)','green',8],\
        #['8c/axon_umbrella_prot8c_glua2_ph','8c: NE2+D1 (acidic)','green',8],\
        #['8d/axon_umbrella_prot8d_glua2_ph','8d: NE2+D1 (acidic)','green',8],\
        #['8e/axon_umbrella_prot8e_glua2_ph','8e: NE2+D1 (acidic)','green',8],\
        #['16a/axon_umbrella_prot16a_glua2_ph','16a: H205A (mutant)','blue',16],\
        #['16b/axon_umbrella_prot16b_glua2_ph','16b: H205A (mutant)','blue',16],\
        #['16c/axon_umbrella_prot16c_glua2_ph','16c: H205A (mutant)','blue',16],\
        #['16d/axon_umbrella_prot16d_glua2_ph','16c: H205A (mutant)','blue',16],\
        #['16e/axon_umbrella_prot16e_glua2_ph','16f: H205A (mutant)','blue',16],\
        #['16f/axon_umbrella_prot16f_glua2_ph','16e: H205A (mutant)','blue',16]\
        ['35/rep1/axon_umbrella_prot35_rep1_glua2_ph','35-1: NE2, default (neutral)','red',35],\
        ['35/rep2/axon_umbrella_prot35_rep2_glua2_ph','35-2: NE2, default (neutral)','red',35],\
        ['35/rep3/axon_umbrella_prot35_rep3_glua2_ph','35-3: NE2, default (neutral)','red',35],\
        ['35/rep4/axon_umbrella_prot35_rep4_glua2_ph','35-4: NE2, default (neutral)','red',35],\
        ['35/rep5/axon_umbrella_prot35_rep5_glua2_ph','35-5: NE2, default (neutral)','red',35],\
        ['36/rep1/axon_umbrella_prot36_rep1_glua2_ph','36-1: NE2, default (neutral)','blue',36],\
        ['36/rep2/axon_umbrella_prot36_rep2_glua2_ph','36-2: NE2, default (neutral)','blue',36],\
        ['36/rep3/axon_umbrella_prot36_rep3_glua2_ph','36-3: NE2, default (neutral)','blue',36],\
        ['36/rep4/axon_umbrella_prot36_rep4_glua2_ph','36-4: NE2, default (neutral)','blue',36],\
        ['36/rep5/axon_umbrella_prot36_rep5_glua2_ph','36-4: NE2, default (neutral)','blue',36]\
        ]

## get indices between dist_min and dist_max
dist_min = 5.0 # nm
dist_converge = 7.5 # nm
path = elements[0][0]
filename = '%s/profile.xvg' % path
dist,E = np.genfromtxt(filename,skip_header=17,skip_footer=4,usecols=[0,1],unpack=True)
idx = np.where((dist >= dist_min) & (dist <= dist_converge))

## min and max dist
dist_min = 5.0 # nm
dist_converge = 7.5 # nm

## loop over elements
folder_previous = 0
list6,list7,list8,list16,list35,list36 = [],[],[],[],[],[]
sum6,sum7,sum8,sum16,sum35,sum36 = 0.0,0.0,0.0,0.0,0.0,0.0
count6,count7,count8,count16,count35,count36 = 0,0,0,0,0,0
dist6,dist7,dist8,dist16,dist35,dist36 = [],[],[],[],[],[]
for e in elements:
    path   = e[0]
    label  = e[1]
    color  = e[2]
    folder = e[3]
    filename = '%s/profile.xvg' % path 
    dist,E = np.genfromtxt(filename,skip_header=17,skip_footer=4,usecols=[0,1],unpack=True)
    
    ## get indices between dist_min and dist_max
    if folder != folder_previous:
        dist,E = np.genfromtxt(filename,skip_header=17,skip_footer=4,usecols=[0,1],unpack=True)
        idx = np.where((dist >= dist_min) & (dist <= dist_converge))
        folder_previous = folder

    ## get energies    
    E_conv = E[idx]
    dist_conv = dist[idx]
    Emin = min(E)
    Emax = np.mean(E_conv[-15:-1])
    Edif = Emax - Emin
    Efin = E_conv-E_conv[-1]
    print('%-34s:  \tEmin,Emax,Edif = %1.1f, %1.1f, %1.1f' % (path,Emin,Emax,Edif))
    
    ## plot curve
    if PLOT_ALL:
        plt.plot(dist_conv,Efin,color=color,alpha=0.3)

    ## add values for later
    if folder == 7:
        sum7 += Efin
        count7 += 1
        list7.append(Efin)
        dist7 = dist_conv
    elif folder == 6:
        sum6 += Efin
        count6 += 1 
        list6.append(Efin)
        dist6 = dist_conv
    elif folder == 8:
        sum8 += Efin
        count8 += 1
        list8.append(Efin)
        dist8 = dist_conv
    elif folder == 16:
        sum16 += Efin
        count16 += 1
        list16.append(Efin)
        dist16 = dist_conv
    elif folder == 35:
        sum35 += Efin
        count35 += 1
        list35.append(Efin)
        dist35 = dist_conv
    elif folder == 36:
        sum36 += Efin
        count36 += 1
        list36.append(Efin)
        dist36 = dist_conv

if count6 > 0:
    matrix6  = np.vstack(list6)
    std6     = matrix6.std(0)
    err6     = std6/np.sqrt(count6)
    mean6    = sum6/count6
    plt.plot(dist6,mean6,color='grey',label='H at D1 (neutral, alternative)')
    plt.fill_between(dist6,mean6-err6,mean6+err6,color='orange',alpha=0.5)

if count7 > 0:
    matrix7  = np.vstack(list7)
    std7     = matrix7.std(0)
    err7     = std7/np.sqrt(count7)
    mean7    = sum7/count7
    plt.plot(dist7,mean7,color='grey',label=r'N$_\mathrm{D1}$ deprotonated (neutral)')
    plt.fill_between(dist7,mean7-err7,mean7+err7,color='grey',alpha=0.5)
    PMF7     = np.mean(mean7[0:25])
    stdlow7  = np.mean(err7[0:25])
    stdhigh7 = np.mean(err7[-14:])
    stdPMF7  = np.sqrt(stdlow7**2 + stdhigh7**2)
    print('PMF deprotonated = %f +/- %f' % (PMF7,stdPMF7))

if count8 > 0:
    matrix8  = np.vstack(list8)
    std8     = matrix8.std(0)
    err8     = std8/np.sqrt(count8)
    mean8    = sum8/count8
    plt.plot(dist8,mean8,color='green',label=r'N$_\mathrm{D1}$ protonated (acidic)')
    plt.fill_between(dist8,mean8-err8,mean8+err8,color='green',alpha=0.5)
    PMF8     = np.mean(mean8[0:7])
    stdlow8  = np.mean(err8[0:7])
    stdhigh8 = np.mean(err8[-14:])
    stdPMF8  = np.sqrt(stdlow8**2 + stdhigh8**2)
    print('PMF protonated = %f +/- %f' % (PMF8,stdPMF8))

if count16 > 0:
    matrix16 = np.vstack(list16)
    std16    = matrix16.std(0)
    err16    = std16/np.sqrt(count16)
    mean16   = sum16/count16
    plt.plot(dist16,mean16,color='blue',label='H212A (mutant)')
    plt.fill_between(dist16,mean16-err16,mean16+err16,color='blue',alpha=0.5)

if count35 > 0:
    matrix35 = np.vstack(list35)
    std35    = matrix35.std(0)
    err35    = std35/np.sqrt(count35)
    mean35   = sum35/count35
    plt.plot(dist35,mean35,color='red',label=r'N$_\mathrm{D1}$ deprotonated (neutral)')
    plt.fill_between(dist35,mean35-err35,mean35+err35,color='red',alpha=0.5)
    PMF35     = np.mean(mean35[0:25])
    stdlow35  = np.mean(err35[0:25])
    stdhigh35 = np.mean(err35[-14:])
    stdPMF35  = np.sqrt(stdlow35**2 + stdhigh35**2)
    print('PMF deprotonated = %f +/- %f' % (PMF35,stdPMF35))

if count36 > 0:
    matrix36 = np.vstack(list36)
    std36    = matrix36.std(0)
    err36    = std36/np.sqrt(count36)
    mean36   = sum36/count36
    plt.plot(dist36,mean36,color='blue',label=r'N$_\mathrm{D1}$ protonated (acidic)')
    plt.fill_between(dist36,mean36-err36,mean36+err36,color='blue',alpha=0.5)
    PMF36     = np.mean(mean36[0:25])
    stdlow36  = np.mean(err36[0:25])
    stdhigh36 = np.mean(err36[-14:])
    stdPMF36  = np.sqrt(stdlow36**2 + stdhigh36**2)
    print('PMF deprotonated = %f +/- %f' % (PMF36,stdPMF36))

dPMF = PMF36-PMF35
std_dPMF = np.sqrt(stdPMF36**2 + stdPMF35**2)
print('Delta PMF = %f +/- %f' % (dPMF,std_dPMF))
xlim=[5.0,7.5]
plt.plot(xlim,[0,0],color='grey',linestyle='--')
plt.xlabel('NTD distance [nm]')
plt.ylabel(r'Free energy, $\Delta G$ [kJ/mol]')
plt.xlim(xlim)
plt.legend(title='His212 protonation state',frameon=False,loc='lower right')
plt.savefig('GluA2_PMF.pdf',format='pdf')
plt.show()

