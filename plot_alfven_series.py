import matplotlib.pyplot as plt
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.animation import FuncAnimation

files_list = [0, 100, 200, 300]


fig = plt.figure()
fig.suptitle(u't_stop = 0.458 s')

ax = fig.add_subplot(1, 2, 1)
ax.set_title(u'V(z)')
ax.set_xlabel(r'$z, см$')
ax.set_ylabel(r'$V, см/с$')

bx = fig.add_subplot(1, 2, 2)
bx.set_title(u'B(z)')
bx.set_xlabel(r'$z, см$')
bx.set_ylabel(r'$B, Гс$')


for file in files_list:
    data1 = np.loadtxt("data_v" + str(file) + ".txt")   
    ax.plot(data1[:,0], data1[:,1], '-', label=str(file))
    bx.plot(data1[:,0], data1[:,2], '-', label=str(file))


ax.legend(loc= 'upper left')

plt.show()