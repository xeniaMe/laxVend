import matplotlib.pyplot as plt
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter
import matplotlib.animation as animation


# настройка размера шрифта на рисунке
fsize = 12
plt.rcParams.update({'axes.titlesize': fsize})
plt.rcParams.update({'axes.labelsize': fsize})
plt.rcParams.update({'font.size': fsize})
plt.rcParams.update({'xtick.labelsize': fsize})
plt.rcParams.update({'ytick.labelsize': fsize})
plt.rcParams.update({'legend.fontsize': fsize-2})

path = "./out/"

file0 = 0 # номер первого файла данных
fileN = 1600 # номер последнего файла данных
dfile = 40 # шаг по файлам - с каким интервалом их считывать
f_list = [] # список файлов для анимации
for i in range(int(fileN/dfile)):
    file_n = (i+1)*dfile
    f_list.append(file_n)
    
# f_list.append(783)
files_list = f_list

# 1 дюйм в сантиметрах
inch2cm = 2.54
# ширина рисунка в см
fig_width_cm = 15
# высота рисунка в см
fig_height_cm = 10

# создание окна рисунка с 2-мя панелями
fig, axs = plt.subplots(1, 1, figsize=(
    fig_width_cm/inch2cm, fig_height_cm/inch2cm), dpi=200)

# панель скорости
# axs.set_title(u'скорость')
axs.set_xlabel(r'$z, км$')
axs.set_ylabel(r'$B, Гс$')
axs.set_ylim(-0.1, 0.1)
axs.grid()

# единица измерения координат
r_unit = 1e5  # 1 километр

ims = []
for file in files_list:
    data1 = np.loadtxt(path + "data_v" + str(file) + ".dat")
    # axs.set_title(label=str(file))

    im = axs.plot(data1[:,0]/r_unit, data1[:,2], '-', color='black', markersize=2.0, linewidth=1, label=str(file))
    ims.append(im)


plt.tight_layout()

ani = animation.ArtistAnimation(fig, ims, interval=500, blit=True,                          repeat_delay=1000)
# ani.save("v.gif", dpi=300, writer=PillowWriter(fps=100))
