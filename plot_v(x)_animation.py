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


file0 = 0 # номер первого файла данных
fileN = 700 # номер последнего файла данных
dfile = 10 # шаг по файлам - с каким интервалом их считывать
f_list = [] # список файлов для анимации
for i in range(int(fileN/dfile)):
    file_n = (i+1)*dfile
    f_list.append(file_n)
    
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
axs.set_ylabel(r'$V, км/с$')

# единица измерения координат
r_unit = 1e5  # 1 километр
# единица измерения скорости
v_unit = 1e5  # 1 километр в секунду


ims = []
for file in files_list:
    data1 = np.loadtxt("data_v" + str(file) + ".txt")
    im = axs.plot(data1[:,0]/r_unit, data1[:,1] / v_unit, '-', label=str(file))
    ims.append(im)


plt.tight_layout()

ani = animation.ArtistAnimation(fig, ims, interval=500, blit=True,                          repeat_delay=1000)
ani.save("v.gif", dpi=300, writer=PillowWriter(fps=100))
