import numpy as np
import math 
from timeit import default_timer as timer

import warnings
warnings.filterwarnings("ignore")

# ----- константы -------------------------------------------------------------
au = 1.5e13        # астрономическая единица, см
G = 6.67e-8        # универсальная грав. постоянная, СГС
Msun = 1.9e33      # масса Слнца, г
Mstar = 1.0 * Msun # масса звезды, в массах Солнца
Rg = 8.31e7        # универсальная газовая постоянная, СГС
mu = 2.3           # молекулярный вес газа
year = 365*24*3600 # 1 год, сек
m_p = 1.67e-24     # масса протона, г

# ----- параметры схемы -------------------------------------------------------
N = 512   # число узлов
Ngs = 1   # яисло фиктивных ячеек с каждой стороны
Ntot = N + 2*Ngs # полное число ячеек
t = 0     # текущее время
c = 0.95  # Число Куранта

L = 10    # длина расчетной области, в единицах шкалы высот диска H
a = 0     # координата левой границы расчетной области
b = L     # координата правой границы расчетной области

eps = 100.0 # параметр искусственной вязкости

# ----- параметры физ. модели -------------------------------------------------
r_cm = 0.052 * au          # радиальное расстояние от звезды, а.е.
z_d = 5.0                  # полутолщина диска, в единицах шкалы высот диска H
n_ISM = 1e9                # концентрация газа в мезвездной среде, см^(-3)
rho_ISM = mu * m_p * n_ISM # плотность МЗС, г/см^3
beta = 80.0                # плазменный параметр в диске
Bz = 1.0                   # безразмерная компонента Bz

# 0 - расчет без вращения, 1 - с вращением
rotation_flag = 1.0 
# 0 - расчет без гравитации, 1 - с гравитацией
gravity_flag = 1.0

# ----- масштабы безразмерных переменных --------------------------------------
rho0 = 1.28e-7   # температура в экватор. плоскости диска, г/см^3
T0 = 1530        # температура внутри диска, К
c_Td = np.sqrt(Rg * T0 / mu) # скорость звука внутри диска, см/с
v0 = c_Td
H = 9e-4*au     # шкала высот диска, см
t0 = H / v0     # шкала измерения времени, сек
p0 = rho0*v0**2 # шкала измерения давления, дин
B0 = np.sqrt(8.0*np.pi*p0/beta) # шкала измерения маг. поля, Гс
g0 = v0**2/H    # шкала измерения ускорения, см/с^2

# ----- коэффициенты уравнений ------------------------------------------------

Ta = 3700                    # температура в атомосфере, К
c_Ta = np.sqrt(Rg * Ta / mu) # скорость звука в атмосфере, см/с

# кеплеровская скорость, см/с
def v_k(r_au):
    return np.sqrt(G * Mstar / (r_au * au))

# кеплеровская скорость в начальной точке, см/с
vk0 = v_k(r_cm / au)

# скорость вращения в случае центробежного равновесия, см/с
def v_phi0(r_cm, z_cm):
    return v_k(r_cm / au) * np.power(1.0 + (z_cm / r_cm)**2, -0.75)

# профиль плотности изотермического гидростатического диска, беразмерная
# здесь z - безразмерная координата (в единицах H)
def rho_disk(z):
    return np.exp(-0.5 * (vk0/c_Td)**2 * (1.0 - (1.0 + (z*H/r_cm)**2.0)**(-0.5)))

rho_surf = rho_disk(z_d)*T0/Ta
# профиль плотности изотермической гидростатической атмосферы диска, беразмерная
# здесь z - безразмерная координата (в единицах H)
def rho_atm(z):
    # плотность на текущей высоте
    rho_at_z = rho_surf * np.exp((vk0/c_Ta)**2 * ((1.0 + (z*H/r_cm)**2.0)**(-0.5) - (1.0 + (z_d*H/r_cm)**2.0)**(-0.5)))
    if rho_at_z > rho_ISM:
        return rho_surf * np.exp((vk0/c_Ta)**2 * ((1.0 + (z*H/r_cm)**2.0)**(-0.5) - (1.0 + (z_d*H/r_cm)**2.0)**(-0.5)))
    else:
        return rho_ISM

# ускорение силы тяжести звезды, см/с^2
def g_z(r_cm, z_cm, M_g):
    return gravity_flag * G*M_g*z_cm/np.power(r_cm**2 + z_cm**2, 1.5)

# альвеновская скорость в начале расчета, см/с
vA = B0 / np.sqrt(4.0 * np.pi * rho0)
# время окончания расчета (безразмерное)
t_stop = 0.2
# начальная максимальная скорость на сетке, см/с
v_max = abs(vA)

print("Model scales:")
print(" rho0   = ", rho0, " g/cm^3")
print(" t0 = ", t0/year, " yr")
print(" v0 = ", v0/1e5, " km/s")
print(" B0     = ", B0, " G")
print(" p0     = ", p0, " dyn")
print(" g0     = ", g0, " dyn")

print("Model characteristics:")
print(" vA     = ", vA / 1e5, " km/s")
print(" v_k    = ", vk0 / 1e5, " km/s")
print(" M_k    = ", vk0 / c_Td)
print(" t_stop = ", t_stop, " [t0]")
print(" v_max  = ", v_max / 1e5, " km/s")

# ----- основные переменные ---------------------------------------------------

# консервативные переменные на шаге t^n
u0_n = np.zeros((Ntot)) # = rho
u1_n = np.zeros((Ntot)) # = rho*v_z
u2_n = np.zeros((Ntot)) # = rho*v_phi
u3_n = np.zeros((Ntot)) # = bphi
u4_n = np.zeros((Ntot)) # = p
# консервативные переменные на шаге t^(n+1/2) (промежуточном)
u0_s = np.zeros((Ntot)) # = rho
u1_s = np.zeros((Ntot)) # = rho*v_z
u2_s = np.zeros((Ntot)) # = rho*v_phi
u3_s = np.zeros((Ntot)) # = bphi
u4_s = np.zeros((Ntot)) # = p
# консервативные переменные на шаге t^(n+1)
u0_n1 = np.zeros((Ntot)) # = rho
u1_n1 = np.zeros((Ntot)) # = rho*v_z
u2_n1 = np.zeros((Ntot)) # = rho*v_phi
u3_n1 = np.zeros((Ntot)) # bphi
u4_n1 = np.zeros((Ntot)) # p

# списки векторов консервативных переменных
u_n  = [u0_n, u1_n, u2_n, u3_n, u4_n]
u_s  = [u0_s, u1_s, u2_s, u3_s, u4_s]
u_n1 = [u0_n1, u1_n1, u2_n1, u3_n1, u4_n1]

# примитивные переменные на шаге t^n
rho_n  = np.zeros((Ntot))  # = rho
vz_n    = np.zeros((Ntot)) # = v_z
vphi_n = np.zeros((Ntot))  # = v_phi
B_n    = np.zeros((Ntot))  # = bphi
p_n    = np.zeros((Ntot))  # = p

# список векторов примитивных переменных
pv  = [rho_n, vphi_n, vz_n, B_n, p_n]
# вектор полной скорости на сетке
vv = np.zeros((N))
# массив узлов сетки
zs = np.linspace(a, b, N)


# ----- построение расчетной сетки --------------------------------------------
dz = (b-a)/ (N-1)
dt = c*dz/(v_max/v0)
print(" dz = ", dz, " [H]")
print(" dt = ", dt, " [t0]")

zs[0] = a
for i in range (1, N):
    zs[i] = zs[i-1] + dz
    
# ----- функци начальных и граничных условий ----------------------------------

# начальное распределение безразмерной температуры
def T_IC(z):
    if z < z_d:
        return 1.0
    else:
        return Ta / T0
    
# начальное распределение безразмерной плотности
def rho_IC(z):
    if z < z_d:
        return rho_disk(z)
    else:
        #return rho_surf#rho_atm(z)
        return rho_atm(z)

# начальное распределние безразмерной скорости v_z
def v_IC(z):
    return 0.0

vphi_s = rotation_flag*v_phi0(r_cm, 1.0 * H) / v0
# начальное распределение безразмерной скорости v_phi
def vphi_IC(z):
    if z < z_d:
        return rotation_flag*v_phi0(r_cm, z * H) / v0
    else:
        return 0.0
    
# начальное распределение безразмерной Bphi    
def B_IC(z):
    return 0.0

# начальное распределение безразмерного давления p
def p_IC(z):
    return rho_IC(z) * T_IC(z)

def p2c(var_n, cv, pv):
    """
    Перевод вектора примитивных переменных в вектор консервативных переменных для переменной с номером var_n = {0, 1, 2, 3, 4}

    Parameters
    ----------
    var_n : целое числое
        номер переменной.
    cv : массив
        вектор консервативных переменных.
    pv : массив
        вектор примитивных переменных.

    Returns
    -------
    None.

    """    
    # pv = {rho, vphi, vz, Bphi, p}
    
    if var_n == 0:
        cv[var_n] = pv[0]
    elif var_n == 1:
        cv[var_n] = pv[0] * pv[1] # rho*v_phi
    elif var_n == 2:
        cv[var_n] = pv[0] * pv[2] # rho*v_z
    elif var_n == 3:
        cv[var_n] = pv[3] # Bphi
    elif var_n == 4:
        cv[var_n] = pv[4] # p
    else:
        None

def c2p(var_n, cv, pv):
    """
     Перевод вектора консервативных переменных в вектор примитивных переменных для переменной с номером var_n = {0, 1, 2, 3, 4}

     Parameters
     
     ----------
     var_n : целое числое
         номер переменной.
     cv : массив
         вектор консервативных переменных.
     pv : массив
         вектор примитивных переменных.

     Returns
     -------
     None.

     """    
    if var_n == 0:
        pv[var_n] = cv[0] # rho
    elif var_n == 1:
        pv[var_n] = cv[1] / cv[0] # v_phi = [rho*v_phi]/rho
    elif var_n == 2:
        pv[var_n] = cv[2] / cv[0]  # v_z= [rho*v_z] / rho
    elif var_n == 3:
        pv[var_n] = cv[3] # Bphi
    elif var_n == 4:
        pv[var_n] = cv[4] # p
    else:
        None

# ----- установка начальных условий -------------------------------------------
def SetIC():
    global t, u_n, pv
    
    t = 0.0
    for i in range(Ngs, Ntot-Ngs): # без учета фиктивных ячеек: i=[Ngs, Ntot-Ngs) = [1, Ntot-2]
        # учитываем, что в массиве узлов zs нет координат фиктивных узлов    
        pv[0][i] = rho_IC(zs[i-Ngs])  # rho_n
        pv[1][i] = vphi_IC(zs[i-Ngs]) # vphi_n
        pv[2][i] = v_IC(zs[i-Ngs])    # vz_n
        pv[3][i] = B_IC(zs[i-Ngs])    # B_n
        pv[4][i] = p_IC(zs[i-Ngs])    # p_n

    # перевод заданных начальных примитивных переменных в начальные консервативные
    for var_n in range(5):
        p2c(var_n, u_n, pv)
            

# ----- функции потоков -------------------------------------------------------

# поток массы в ячейке с индексом i
def F0(u, i):
    # = rho*v_z
    return u[2][i]

# поток импульса rho*v_phi в ячейке с индексом i
def F1(u, i):
    # = rho*v_z*v_ph -  (2/beta)*Bz*Bphi
    return u[1][i] * u[2][i] / u[0][i] - 2 * u[3][i] * Bz / beta

# поток импульса rho*v_z в ячейке с индексом i
def F2(u, i):
    # = rho*v_z^2 + p + B_phi^2/beta
    return (u[2][i])**2 / u[0][i] + u[4][i] + (u[3][i])**2/beta

# поток Bphi в ячейке с индексом i
def F3(u, i):
    # = v_z*B_phi - v_phi*Bz
    return u[2][i] * u[3][i] / u[0][i] - u[1][i] * Bz / u[0][i]

# поток давления в ячейке с индексом i
def F4(u, i):
    # = (rho*v_z) * p / rho
    return u[2][i] * u[4][i] / u[0][i]

# поток величины с индексом var_n в ячейке с индексом i
def F(u, i, var_n):
    if var_n == 0:
        return F0(u, i)
    elif var_n == 1:
        return F1(u, i)
    elif var_n == 2:
        return F2(u, i)
    elif var_n == 3:
        return F3(u, i)
    elif var_n == 4:
        return F4(u, i)
    else:
        None

# вектор-функция источников (правые части уравнений)
def Source(u, i, var_n):
    if var_n == 0:
        return 0.0
    elif var_n == 1:
        return 0.0
    elif var_n == 2:
        # безразмерная сила тяжести
        if (i == Ntot-Ngs-1): # этой точки нет в массиве zs
            return -g_z(r_cm, (zs[i-1] + dz)*H, Mstar)/g0 * u[0][i]
        else:
            return -g_z(r_cm, zs[i]*H, Mstar)/g0 * u[0][i]
        
    elif var_n == 3:
        return 0.0
    elif var_n == 4:
        return 0.0
    else:
        None
        
# установка граничных условий
def SetBC():
    global t, u_n1

    # левая граница
    rho_L  = 1.0 # 
    vphi_L = rotation_flag*vk0 / v0 #
    vz_L   = 0.0 #
    B_L    = 0.0 #
    p_L    = 1.0 # 
    
    # это значения на левой границе (i = Ngs = 1)
    u_n1[0][Ngs] = rho_L
    u_n1[1][Ngs] = rho_L * vphi_L
    u_n1[2][Ngs] = rho_L * vz_L
    u_n1[3][Ngs] = B_L
    u_n1[4][Ngs] = p_L
  

    # условия свободного втекания, посредством фиктивных ячеек
    u_n[0][Ntot-1] = u_n[0][Ntot-Ngs-2]
    u_n[1][Ntot-1] = u_n[1][Ntot-Ngs-2]
    u_n[2][Ntot-1] = u_n[2][Ntot-Ngs-2]
    u_n[3][Ntot-1] = u_n[3][Ntot-Ngs-2]
    u_n[4][Ntot-1] = u_n[4][Ntot-Ngs-2]

    # u_n1[0][Ntot-Ngs-1] = rho_R
    # u_n1[1][Ntot-Ngs-1] = rho_R * vphi_R
    # u_n1[2][Ntot-Ngs-1] = rho_R * vz_R
    # u_n1[3][Ntot-Ngs-1] = B_R
    

# вычисление шага по времени
def UpdateTimeStep():
    global dt, v_max, v, v_n
    # флаг указывает, успешно ли выполнен расчет шага
    success = True
    # сообщение об успешности (неуспешности) расчета
    message = "Time step is successfully determined"

    for i in range(Ngs, Ntot-Ngs): # цикл без учета фиктивных ячеек
        # квадрат безразмерной альв. скорости
        vAvA = (2.0/beta) * (Bz**2 + pv[3][i]**2) / pv[0][i]
        # полная скорость, 1.0 - это б.р. скорость звука
        #vv[i-Ngs] = abs(pv[2][i]) + abs(pv[1][i]) + np.sqrt(1.0 + vAvA)
        vv[i-Ngs] = max(pv[2][i], abs(pv[1][i]), np.sqrt(1.0 + vAvA))
        
    v_max = max(vv)
    dt = c*dz/v_max

    if (v_max > 10000.0*vA):
        success = False
        message = ("Velocity became unphysically large: v_max = %.3e" % (v_max/1e5))
    
    return [success, message]


# шаг по времени
def Step():
    global dt, u_n, u_s, u_n1, dz
    
    # Метод Л-В, этап предиктора
    for var_n in range(5):
        for i in range(Ngs, Ntot-Ngs):
            u_s[var_n][i] =  0.5 * (u_n[var_n][i+1] + u_n[var_n][i]) - 0.5*(dt/dz)*(F(u_n, i+1, var_n) - F(u_n, i, var_n)) + Source(u_n, i, var_n)*dt
            
    # Метод Л-В, этап корректора
    for var_n in range(5):
        for i in range (Ngs + 1, Ntot-Ngs):
            u_n1[var_n][i] = u_n[var_n][i] - (dt/dz) * (F(u_s, i, var_n) - F(u_s, i-1, var_n)) + 0.5 * dt * (Source(u_s, i, var_n) + Source(u_n, i, var_n)) + eps * dt * (u_n[var_n][i+1] - 2*u_n[var_n][i] + u_n[var_n][i-1])
            
    # построить на основе решения вектор примитивных переменных
    for var_n in range(5):
        c2p(var_n, u_n1, pv)

# обновление НУ
def UpdateIC():
    global t, u_n, u_n1
    for var_n in range(5):
        # перебируем значение внутри сетки, без учета фиктивных ячеек
        u_n[var_n][Ngs:Ntot-Ngs:] = u_n1[var_n][Ngs:Ntot-Ngs:]

# сохранение данных в файл
def SaveData(n):
    data = []
    data.append(zs)
    
    if (n == 0): # начальное условие
        for var_n in range(5):
            c2p(var_n, u_n, pv)
            # сохраняем данные без учета фиктивных ячеек
            data.append(pv[var_n][Ngs:Ntot-Ngs:])

    else:
        for var_n in range(5):
            c2p(var_n, u_n1, pv)
            # сохраняем данные без учета фиктивных ячеек
            data.append(pv[var_n][Ngs:Ntot-Ngs:])
 
    # столбцы 1, 2, 3, 4, 5, 6 - (безразмерные) массивы zs, rho, vphi, vz, Bphi, p соответственно
    np.savetxt("./out/data" + str(n) + ".dat", np.array(data).transpose(), fmt=('%.3e', '%.3e', '%.3e', '%.3e', '%.3e', '%.3e'))


print("Press Enter to start a simulation")
input()

n = 0 # номер шага по времени
SetIC()
# Сохранение НУ
SaveData(n)

n = n + 1
while t <= t_stop:
    start = timer()
    [contin, message] = UpdateTimeStep()  
    # выход из цикла и сохранение результатов, если что-то пошло не так (TO DO)
    if(contin == False):
        print("Exiting because", message)
        SaveData(n)
        break

    SetBC()
    
    Step()
    # вывод n, t, dt на экран и сохранение результатов каждые ... шагов
    if ((n % 500) == 0):
        print(  "step No ", n)
        print(  "t     = %.3e [t0] = %.1e [s]" % (t, t*t0))
        print(  "dt    = %.3e [t0] = %.1e [s]" % (dt, dt*t0))
        print(  "v_max = %.2e  [km/s]" % (v_max/1e5))
        SaveData(n)
    UpdateIC()  
    t += dt
    n += 1

end = timer()
print(u'Время расчета: %.3e сек' % (end-start))
print("Total number of steps = ", n)
SaveData(n)