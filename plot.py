import matplotlib.pyplot as plt
import numpy as np
from matplotlib import pyplot as plt
import numpy as np
data = np.loadtxt("data.txt")
fig = plt.figure()

ax = fig.add_subplot(1, 1, 1)
#ax.set_xlim(0, 1.0)
#ax.set_ylim(0, 40.0)
ax.set_xlabel(r'$x$')
ax.set_ylabel(r'$u$')
ax.set_title(u'Данные из файла')
ax.plot(data[:,0], data[:,1], 'ko-', label=r'$d_1$')
ax.legend(loc='best')
plt.show()