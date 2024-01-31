import matplotlib.pyplot as plt
import numpy as np
from matplotlib import pyplot as plt
import numpy as np
data1 = np.loadtxt("data_c_1_05.txt")
data2 = np.loadtxt("data_c_0_9.txt")
data3 = np.loadtxt("data_c_0_6.txt")
data4 = np.loadtxt("data_c_0_3.txt")

data_n_20 = np.loadtxt("C:data_n_20.txt")
data_n_80 = np.loadtxt("C:data_n_80.txt")
data_n_160 = np.loadtxt("data_n_160.txt")
data_n_320 = np.loadtxt("data_n_320.txt")
fig = plt.figure()

ax = fig.add_subplot(1, 2, 1)
ax.set_xlabel(r'$x$')
ax.set_ylabel(r'$u$')
ax.set_title(u'Зависимость от числа Куранта при N=40')
ax.plot(data1[:,0], data1[:,1], 'k-', color = 'red', label=r'$c = 1.05$')
ax.plot(data2[:,0], data2[:,1], 'k-', color = 'green', label=r'$c = 0.9$')
ax.plot(data3[:,0], data3[:,1], 'k-', color = 'blue',label=r'$c = 0.6$')
ax.plot(data4[:,0], data4[:,1], 'k-', label=r'$c = 0.3$')
ax.legend(loc='best')

bx = fig.add_subplot(1, 2, 2)
bx.set_xlabel(r'$x$')
bx.set_ylabel(r'$u$')
bx.set_title(u'Зависимость от числа узлов при С = 0.9')
bx.plot(data_n_20[:,0], data_n_20[:,1], 'k-', color = 'red', label=r'$N = 20$')
bx.plot(data2[:,0], data2[:,1], 'k-', color = 'green', label=r'$N = 40$')
bx.plot(data_n_80[:,0], data_n_80[:,1], 'k-', color = 'blue',label=r'$N = 80$')
bx.plot(data_n_160[:,0], data_n_160[:,1], 'k-', label=r'$N = 160$')
bx.plot(data_n_320[:,0], data_n_320[:,1], 'k-', color = 'pink',label=r'$N = 320$')
bx.legend(loc='best')
plt.show()