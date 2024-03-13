import matplotlib.pyplot as plt
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.animation import FuncAnimation

files_list = [0, 100, 200, 300]


fig = plt.figure()
fig.suptitle(u't_stop = 4.5 s')

ax = fig.add_subplot(1, 2, 1)
ax.set_title(u'v(x)')
ax.set_xlabel(r'$x$')
ax.set_ylabel(r'$v$')

bx = fig.add_subplot(1, 2, 2)
bx.set_title(u'b(x)')
bx.set_xlabel(r'$x$')
bx.set_ylabel(r'$b$')


for file in files_list:
    data1 = np.loadtxt("data_v" + str(file) + ".txt")    
    data2 = np.loadtxt("data_b" + str(file) + ".txt")  
    ax.plot(data1[:,0], data1[:,1], '-', label=str(file))
    bx.plot(data2[:,0], data2[:,1], '-', label=str(file))


ax.legend(loc= 'upper left')

plt.show()