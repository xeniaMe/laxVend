import matplotlib.pyplot as plt
import numpy as np
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

lines = []
for _ in files_list:
    line, = ax.plot([], [], '-')
    lines.append(line)

    line, = bx.plot([], [], '-')
    lines.append(line)


def init():
    for line in lines:
        line.set_data([], [])
    return lines

def animate(i):
    for file in files_list:
        data1 = np.loadtxt("data_v" + str(file) + ".txt")
        data2 = np.loadtxt("data_b" + str(file) + ".txt")
        lines[i].set_data(data1[:,0], data1[:,1])
        lines[i].set_data(data2[:,0], data2[:,1])
    return line

# Создаем анимацию
ani = FuncAnimation(fig, animate, init_func=init, frames=200, interval=200, blit=True)

plt.show()