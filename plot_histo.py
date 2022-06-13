import numpy as np
import matplotlib.pyplot as plt

folder = 6
file = '%d/axon_umbrella_prot%d_glua2_ph/histo.xvg' % (folder,folder)

d = np.genfromtxt(file,skip_header=17,usecols=[0],unpack=True)

for i in range(100):
    try:
        h = np.genfromtxt(file,skip_header=17,usecols=[i+1],unpack=True)
        plt.plot(d,h)
    except:
        pass

plt.xlabel('distance [nm]')
plt.savefig('histo%d.png' % folder)
plt.show()
