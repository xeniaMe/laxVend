import matplotlib.pyplot as plt
import numpy as np
from matplotlib import pyplot as plt
import numpy as np

data = np.loadtxt("data.txt")

fig = plt.figure()

ax = fig.add_subplot(1, 1, 1)
ax.set_xlabel(r'$x$')
ax.set_ylabel(r'$u$')
ax.plot(data[:,0], data[:,1], '-', color = 'red', label=r'$c = 1.05$')
ax.legend(loc='best')

plt.show()