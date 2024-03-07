import matplotlib.pyplot as plt
import numpy as np
from matplotlib import pyplot as plt
import numpy as np

files_list = [0, 100, 200, 300]


fig = plt.figure()

ax = fig.add_subplot(1, 1, 1)
ax.set_title(u'v(x)')
ax.set_xlabel(r'$x$')
ax.set_ylabel(r'$v$')

for file in files_list:
    data1 = np.loadtxt("data" + str(file) + ".txt")    
    ax.plot(data1[:,0], data1[:,1], '-', label=str(file))

ax.legend(loc= 'upper left')

plt.show()