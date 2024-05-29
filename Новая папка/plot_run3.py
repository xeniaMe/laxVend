import matplotlib.pyplot as plt
import numpy as np
from matplotlib import pyplot as plt
import numpy as np

data1 = np.loadtxt("data10.txt")
data2 = np.loadtxt("data11.txt")
data3 = np.loadtxt("data12.txt")
fig = plt.figure()

ax = fig.add_subplot(1, 1, 1)
ax.set_title(u'Решение уравнения переноса для разных моментов времени при N = 80 C=0.9')
ax.set_xlabel(r'$x$')
ax.set_ylabel(r'$u$')
ax.plot(data1[:,0], data1[:,1], '-', color = 'red', label=r'$t = 0.05$')
ax.plot(data2[:,0], data2[:,1], '-', color = 'blue', label=r'$t = 0.25$')
ax.plot(data3[:,0], data3[:,1], '-', color = 'orange', label=r'$t = 0.4$')
ax.legend(loc='best')



plt.show()