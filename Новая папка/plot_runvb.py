import matplotlib.pyplot as plt
import numpy as np
from matplotlib import pyplot as plt
import numpy as np

data1 = np.loadtxt("data305.txt")
data2 = np.loadtxt("data302.txt")
data3 = np.loadtxt("data303.txt")
data4 = np.loadtxt("data4.txt")
fig = plt.figure()

ax = fig.add_subplot(1, 1, 1)
ax.set_title(u'v(x)')
ax.set_xlabel(r'$x$')
ax.set_ylabel(r'$v$')
ax.plot(data1[:,0], data1[:,1], '-', color = 'red', label=r'$ a_v = 0.7, c = 0.95, N = 320 $')
#ax.plot(data2[:,0], data2[:,1], '-', color = 'blue', label=r'$ a_v = 0.7, C = 0.3, N = 320, t_stop = 0.4$')
# ax.plot(data3[:,0], data3[:,1], '-', color = 'orange', label=r'$ a_v = 0.7, C = 1.01, N = 320, t_stop = 0.4$')
#ax.plot(data4[:,0], data4[:,1], '-', color = 'black', label=r'$c = 1.05$')
ax.legend(loc= 'upper left')

plt.show()