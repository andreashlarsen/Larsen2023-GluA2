import numpy as np
import matplotlib.pyplot as plt

t,d,logw = np.genfromtxt('WEIGHTS',skip_header=1,usecols=[0,1,2], unpack=True)
logw -= 740 # done to avoid inf when taking exp of logw
w = np.exp(logw)

with open('weight','w') as f:
    #f.write('weight\n')
    for i in range(len(w)):
        f.write('%e\n' % w[i])

with open('distance','w') as f:
    for i in range(len(d)):
        f.write('%f\n' % d[i])

p1=plt.subplot(311)
p2=plt.subplot(312)
p3=plt.subplot(313)

p1.plot(t,d,color='black')
p1.set_ylabel('distance')

p2.plot(t,logw,color='red')
p2.set_ylabel('log(weight)')

p3.plot(t,w,color='blue')
p3.set_ylabel('weight')
p3.set_xlabel('time')

plt.tight_layout()

plt.show()
