import numpy as np
import matplotlib.pyplot as plt

PLOT_metaD_BIAS = 0
PLOT_WALL_BIAS = 0

## import and convert data
time,dist,bias,wall = np.genfromtxt('COLVAR',comments='#',usecols=[0,1,2,3],unpack=True)
time = time/1000 # ns -> us

## get initial dist
dist_ini = 6 #dist[0]

## define dist for bulgy conformation
dist_bulgy = 14

## plot
rows = 1
if PLOT_WALL_BIAS:
    rows += 1
if PLOT_metaD_BIAS:
    rows += 1

p1 = plt.subplot(rows,1,1)
p1.plot(time,dist,color='black')
p1.plot([time[0],time[-1]],[dist_ini,dist_ini],linestyle='--',color='grey')
p1.plot([time[0],time[-1]],[dist_bulgy,dist_bulgy],linestyle='--',color='grey')
p1.set_ylabel('NTD distance [nm]')
p1.set_xlim(time[0],time[-1])

if PLOT_metaD_BIAS:
    p2 = plt.subplot(rows,1,2)
    p2.plot(time,bias,color='black')
    p2.set_ylabel('MetaD bias [kJ/mol]')
    p2.set_xlim(time[0],time[-1])
    if PLOT_WALL_BIAS:
        p3 = plt.subplot(rows,1,rows)
        p3.plot(time,wall,color='black')
        p3.set_ylabel('Wall bias [kJ/mol]')
        p3.set_xlabel(r'Time [$\mu$s]')
        p3.set_xlim(time[0],time[-1])
    else: 
        p2.set_xlabel(r'Time [$\mu$s]')
        p2.set_xlim(time[0],time[-1])
elif PLOT_WALL_BIAS:
    p2 = plt.subplot(rows,1,rows)
    p2.plot(time,wall,color='black')
    p2.set_ylabel('Wall bias [kJ/mol]')
    p2.set_xlabel(r'Time [$\mu$s]')
    p2.set_xlim(time[0],time[-1])
else:
    p1.set_xlabel(r'Time [$\mu$s]')
    p1.set_xlim(time[0],time[-1])

plt.tight_layout()
plt.savefig('COLVAR.pdf',format='pdf')
plt.show()

