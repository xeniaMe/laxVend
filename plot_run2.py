import matplotlib.pyplot as plt
import numpy as np
from matplotlib import pyplot as plt
import numpy as np

data1 = np.loadtxt("data5.txt")
data2 = np.loadtxt("data6.txt")
data3 = np.loadtxt("data7.txt")
data4 = np.loadtxt("data8.txt")
data5 = np.loadtxt("data9.txt")
fig = plt.figure()

ax = fig.add_subplot(1, 1, 1)
ax.set_title(u'Зависимость от числа узлов при С = 0.9')
ax.set_xlabel(r'$x$')
ax.set_ylabel(r'$u$')
ax.plot(data1[:,0], data1[:,1], '-', color = 'red', label=r'$N = 20$')
ax.plot(data2[:,0], data2[:,1], '-', color = 'blue', label=r'$N = 40$')
ax.plot(data3[:,0], data3[:,1], '-', color = 'orange', label=r'$N = 80$')
ax.plot(data4[:,0], data4[:,1], '-', color = 'black', label=r'$N = 160$')
ax.plot(data5[:,0], data5[:,1], '-', color = 'green', label=r'$N = 320$')
ax.legend(loc='best')



plt.show()