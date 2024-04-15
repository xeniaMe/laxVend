import matplotlib.pyplot as plt
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.animation import FuncAnimation


# настройка размера шрифта на рисунке
fsize=12
plt.rcParams.update({'axes.titlesize':fsize})
plt.rcParams.update({'axes.labelsize':fsize})
plt.rcParams.update({'font.size':fsize})
plt.rcParams.update({'xtick.labelsize':fsize})
plt.rcParams.update({'ytick.labelsize':fsize})
plt.rcParams.update({'legend.fontsize':fsize-2})


#files_list = [0, 100, 200, 300]
files_list = [0, 100]


# 1 дюйм в сантиметрах
inch2cm = 2.54
# ширина рисунка в см
fig_width_cm = 15
# высота рисунка в см
fig_height_cm = 10

# создание окна рисунка с 2-мя панелями
fig, axs = plt.subplots(1, 2, figsize=(fig_width_cm/inch2cm, fig_height_cm/inch2cm), dpi=200)

# панель скорости
axs[0].set_title(u'скорость')
axs[0].set_xlabel(r'$z, км$')
axs[0].set_ylabel(r'$V, км/с$')

# панель магнитного олполя
axs[1].set_title(u'магнитное поле')
axs[1].set_xlabel(r'$z, км$')
axs[1].set_ylabel(r'$B, Гс$')

# единица измерения координат
r_unit = 1e5 # 1 километр
# единица измерения скорости
v_unit = 1e5 # 1 километр в секунду

for file in files_list:
    data1 = np.loadtxt("data_v" + str(file) + ".txt")   
    axs[0].plot(data1[:,0]/r_unit, data1[:,1] / v_unit, '-', label=str(file))
    axs[1].plot(data1[:,0]/r_unit, data1[:,2] / v_unit, '-', label=str(file))


axs[0].legend(loc= 'upper left')
plt.tight_layout()

fig.show()