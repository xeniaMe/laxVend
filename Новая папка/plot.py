import matplotlib.pyplot as plt
import numpy as np
from matplotlib import pyplot as plt
import numpy as np
data1 = np.loadtxt("C:\Рабочий стол\work1\data_c_1_05.txt")
data2 = np.loadtxt("C:\Рабочий стол\work1\data_c_0_9.txt")
data3 = np.loadtxt("C:\Рабочий стол\work1\data_c_0_6.txt")
data4 = np.loadtxt("C:\Рабочий стол\work1\data_c_0_3.txt")

data_n_20 = np.loadtxt("C:\Рабочий стол\work1\data_n_20.txt")
data_n_80 = np.loadtxt("C:\Рабочий стол\work1\data_n_80.txt")
data_n_160 = np.loadtxt("C:\Рабочий стол\work1\data_n_160.txt")
data_n_320 = np.loadtxt("C:\Рабочий стол\work1\data_n_320.txt")

data0_05_t = np.loadtxt("C:\Рабочий стол\work1\data0_05_t.txt")
data0_25_t = np.loadtxt("C:\Рабочий стол\work1\data0_25_t.txt")
data0_4_t = np.loadtxt("C:\Рабочий стол\work1\data0_4_t.txt")

fig = plt.figure()

ax = fig.add_subplot(1, 2, 1)
ax.set_xlabel(r'$x$')
ax.set_ylabel(r'$u$')
ax.set_title(u'Зависимость от числа Куранта при N=40')
ax.plot(data1[:,0], data1[:,1], '-', color = 'red', label=r'$c = 1.05$')
ax.plot(data2[:,0], data2[:,1], '-', color = 'green', label=r'$c = 0.9$')
ax.plot(data3[:,0], data3[:,1], '-', color = 'blue',label=r'$c = 0.6$')
ax.plot(data4[:,0], data4[:,1], '-', label=r'$c = 0.3$')
ax.legend(loc='best')

bx = fig.add_subplot(1, 2, 2)
bx.set_xlabel(r'$x$')
bx.set_ylabel(r'$u$')
bx.set_title(u'Зависимость от числа узлов при С = 0.9')
bx.plot(data_n_20[:,0], data_n_20[:,1], '-', color = 'red', label=r'$N = 20$')
bx.plot(data2[:,0], data2[:,1], '-', color = 'green', label=r'$N = 40$')
bx.plot(data_n_80[:,0], data_n_80[:,1], '-', color = 'blue',label=r'$N = 80$')
bx.plot(data_n_160[:,0], data_n_160[:,1], '-', label=r'$N = 160$')
bx.plot(data_n_320[:,0], data_n_320[:,1], '-', color = 'pink',label=r'$N = 320$')
bx.legend(loc='best')

fig2 = plt.figure()
сx = fig2.add_subplot(1, 1, 1)
сx.set_xlabel(r'$x$')
сx.set_ylabel(r'$u$')
сx.set_title(u'Решение уравнения переноса для разных моментов времени при N = 80 C=0.9')
сx.plot(data0_05_t[:,0], data0_05_t[:,1], '-', color = 'r', label=r'$t = 0.05$')
сx.plot(data0_25_t[:,0], data0_25_t[:,1], '-', color = 'g', label=r'$t = 0.25$')
сx.plot(data0_4_t[:,0], data0_4_t[:,1], '-', color = 'b',label=r'$t = 0.4$')
сx.legend(loc='best')
plt.show()