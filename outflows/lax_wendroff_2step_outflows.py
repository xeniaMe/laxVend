import numpy as np
import math 

import matplotlib.pyplot as plt 
import warnings
warnings.filterwarnings("ignore")

# ----- константы -------------------------------------------------------------
au = 1.5e13
G = 6.67e-8
Mstar = 1.99e33
Rg = 8.31e7
mu = 2.3
year = 365*24*3600

# ----- параметры схемы -------------------------------------------------------
N = 512 # число узлов
Ngs = 1 # яисло фиктивных ячеек с каждой стороны
Ntot = N + 2*Ngs # полное число ячеек
t = 0 #текущее время
c = 0.95 #Kurant number

L = 5   # длина расчетной области , H
a = 0   # координата левой границы расчетной области
b = L   # координата правой границы расчетной области

# ----- параметры физ. модели -------------------------------------------------
r_cm = 1 * au # радиальное расстояние от звезды
z_d = 1.0     # полутолщина диска
rho_ISM = np.exp(-0.5*z_d**2) # безразмерная плотность МЗС, в ед. rho0
beta = 1.0    # плазменный параметр
Bz = 1.0      # безразмерная компонента Bz (для нее уравнения не решаются)
eps = 0.2

# 0 - расчет без вращения, 1 - с вращением
rotation_flag = 1.0 
# 0 - расчет без гравитации, 1 - с гравитацией
gravity_flag = 1.0

# ----- масштабы безразмерных переменных --------------------------------------
rho0 = 1e-13 # температура в экватор. плоскости, г/см^3
T0 = 500     # температура, К
c_T = np.sqrt(Rg * T0 / mu) # скорость звука, см/с
v0 = c_T
H = 0.05 * r_cm # шкала высот диска
t0 = H / v0     # шкала измерения времени
beta= 1.0        # плазменный параметр
p0 = rho0*v0**2 # шкала измерения давления
B0 = np.sqrt(8.0*np.pi*p0/beta) # шкала измерения маг. поля
g0 = v0**2/H    # шкала измерения ускорения

# ----- коэффициенты уравнений ------------------------------------------------

# кеплеровская скорость
def v_k(r_au):
    return np.sqrt(G * Mstar / (r_au * au))

vk0 = v_k(r_cm / au)

# скорость вращения в случае центробежного равновесия
def v_phi0(r_cm, z_cm):
    return v_k(r_cm / au) * np.power(1.0 + (z_cm / r_cm)**2, -0.75)

# профиль плотности изотермической гидростатической атмосферы
def rho_hs(z):
    # return np.exp(-0.5 * z**2)
    return np.exp(-0.5 * (vk0/c_T)**2 * (1.0 - (1.0 + (z*H/r_cm)**2.0)**(-0.5)))

# завивимость безразмерного давления от безразмерной плотности
def p(rho):
    return rho

# ускорение силы тяжести звезды, см/с^2
def g_z(r_cm, z_cm, M_g):
    return gravity_flag * G*M_g*z_cm/np.power(r_cm**2 + z_cm**2, 1.5)

# альвеновская скорость в начале расчета
vA = B0 / np.sqrt(4.0 * np.pi * rho0)
# время окончания расчета (безразмерное)
t_stop = 0.25
# начальная максимальная скорость на сетке
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
print(" M_k    = ", vk0 / c_T)
print(" t_stop = ", t_stop, " [t0]")
print(" v_max  = ", v_max / 1e5, " km/s")

# ----- основные переменные ---------------------------------------------------

# консервативные переменные на шаге t^n
u0_n = np.zeros((Ntot)) # = rho
u1_n = np.zeros((Ntot)) # = rho*v_z
u2_n = np.zeros((Ntot)) # = rho*v_phi
u3_n = np.zeros((Ntot)) # = bphi
# консервативные переменные на шаге t^(n+1/2) (промежуточном)
u0_s = np.zeros((Ntot)) # = rho
u1_s = np.zeros((Ntot)) # = rho*v_z
u2_s = np.zeros((Ntot)) # = rho*v_phi
u3_s = np.zeros((Ntot)) # = bphi
# консервативные переменные на шаге t^(n+1)
u0_n1 = np.zeros((Ntot)) # = rho
u1_n1 = np.zeros((Ntot)) # = rho*v_z
u2_n1 = np.zeros((Ntot)) # = rho*v_phi
u3_n1 = np.zeros((Ntot)) # bphi

# списки векторов консервативных переменных
u_n  = [u0_n, u1_n, u2_n, u3_n]
u_s  = [u0_s, u1_s, u2_s, u3_s]
u_n1 = [u0_n1, u1_n1, u2_n1, u3_n1]

# примитивные переменные на шаге t^n
rho_n  = np.zeros((Ntot))  # = rho
vz_n    = np.zeros((Ntot)) # = v_z
vphi_n = np.zeros((Ntot))  # = v_phi
B_n    = np.zeros((Ntot))  # = bphi

# список векторов примитивных переменных
pv  = [rho_n, vphi_n, vz_n, B_n]
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

# начальное распределение плотности
def rho_IC(z):
    if z < z_d:
        return rho_hs(z)#np.exp(-0.5 * z**2)
    else:
        return rho_hs(z_d)#np.exp(-0.5)

# начальное распределние скорости v_z
def v_IC(z):
    return 0.0

vphi_s = rotation_flag*v_phi0(r_cm, 1.0 * H) / v0
# начальное распределение скорости v_phi
def vphi_IC(z):
    if z < 1.0:
        return rotation_flag*v_phi0(r_cm, z * H) / v0
    elif z < 1.3: # плавный переход к атмосфере диска - чтобы избежать распада разрыва
        return vphi_s + (z - 1.0) * (0 - vphi_s) / (1.2 - 1.0)
    else:
        return 0.0
    
# начальное распределение Bphi    
def B_IC(z):
    return 0.0

def p2c(var_n, cv, pv):
    """
    Перевод вектора примитивных переменных в вектор консервативных переменных для переменной с номером var_n = {0, 1, 2, 3}

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
    # pv = {rho, vphi, vz, Bphi}
    
    if var_n == 0:
        cv[var_n] = pv[0]
    elif var_n == 1:
        cv[var_n] = pv[0] * pv[1] # rho*v_phi
    elif var_n == 2:
        cv[var_n] = pv[0] * pv[2] # rho*v_z
    elif var_n == 3:
        cv[var_n] = pv[3] # Bphi
    else:
        None

def c2p(var_n, cv, pv):
    """
     Перевод вектора консервативных переменных в вектор примитивных переменных для переменной с номером var_n = {0, 1, 2, 3}

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

    # перевод заданных начальных примитивных переменных в начальные консервативные
    for var_n in range(4):
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
    return (u[2][i])**2 / u[0][i] + p(u[0][i]) + (u[3][i])**2/beta

# поток Bphi в ячейке с индексом i
def F3(u, i):
    # = v_z*B_phi - v_phi*Bz
    return u[2][i] * u[3][i] / u[0][i] - u[1][i] * Bz / u[0][i]

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
    else:
        None
        
# установка граничных условий
def SetBC():
    global t, u_n1

    # левая граница
    rho_L  = 1.0 # 
    vphi_L = rotation_flag*vk0 / v0 #
    vz_L    = 0.0 #
    B_L    = 0.0 #
    
    # это значения на левой границе (i = Ngs = 1)
    u_n1[0][Ngs] = rho_L
    u_n1[1][Ngs] = rho_L * vphi_L
    u_n1[2][Ngs] = rho_L * vz_L
    u_n1[3][Ngs] = B_L
  

    # условия свободного втекания, посредством фиктивных ячеек
    u_n[0][Ntot-1] = u_n[0][Ntot-Ngs-2]
    u_n[1][Ntot-1] = u_n[1][Ntot-Ngs-2]
    u_n[2][Ntot-1] = u_n[2][Ntot-Ngs-2]
    u_n[3][Ntot-1] = u_n[3][Ntot-Ngs-2]

    # u_n1[0][Ntot-Ngs-1] = rho_R
    # u_n1[1][Ntot-Ngs-1] = rho_R * vphi_R
    # u_n1[2][Ntot-Ngs-1] = rho_R * vz_R
    # u_n1[3][Ntot-Ngs-1] = B_R
    

# вычисление шага по времени
def UpdateTimeStep():
    global dt, v_max, v, v_n

    for i in range(Ngs, Ntot-Ngs): # цикл без учета фиктивных ячеек
        # квадрат безразмерной альв. скорости
        vAvA = (2.0/beta) * (Bz**2 + pv[3][i]**2) / pv[0][i]
        # полная скорость, 1.0 - это б.р. скорость звука
        vv[i-Ngs] = abs(pv[2][i]) + abs(pv[1][i]) + np.sqrt(1.0 + vAvA)
        
    v_max = max(vv)
    dt = c*dz/v_max
    
    return True


# шаг по времени
def Step():
    global dt, u_n, u_s, u_n1, dz
    
    # Метод Л-В, этап предиктора
    for var_n in range(4):
        for i in range(Ngs, Ntot-Ngs):
            u_s[var_n][i] =  0.5 * (u_n[var_n][i+1] + u_n[var_n][i]) - 0.5*(dt/dz)*(F(u_n, i+1, var_n) - F(u_n, i, var_n)) + Source(u_n, i, var_n)*dt
            
    # Метод Л-В, этап корректора
    for var_n in range(4):
        for i in range (Ngs + 1, Ntot-Ngs):
            u_n1[var_n][i] = u_n[var_n][i] - (dt/dz) * (F(u_s, i, var_n) - F(u_s, i-1, var_n)) + 0.5 * dt * (Source(u_s, i, var_n) + Source(u_n, i, var_n)) + eps*(u_n[var_n][i+1]+2*u_n[var_n][i]+u_n[var_n][i-1])
            
    # построить на основе решения вектор примитивных переменных
    for var_n in range(4):
        c2p(var_n, u_n1, pv)

# обновление НУ
def UpdateIC():
    global t, u_n, u_n1
    for var_n in range(4):
        # перебируем значение внутри сетки, без учета фиктивных ячеек
        u_n[var_n][Ngs:Ntot-Ngs:] = u_n1[var_n][Ngs:Ntot-Ngs:]

# сохранение данных в файл
def SaveData(n):
    data = []
    data.append(zs)
    
    if (n == 0): # начальное условие
        for var_n in range(4):
            c2p(var_n, u_n, pv)
            # сохраняем данные без учета фиктивных ячеек
            data.append(pv[var_n][Ngs:Ntot-Ngs:])

    else:
        for var_n in range(4):
            c2p(var_n, u_n1, pv)
            # сохраняем данные без учета фиктивных ячеек
            data.append(pv[var_n][Ngs:Ntot-Ngs:])
 
    # столбцы 1, 2, 3, 4, 5 - (безразмерные) массивы zs, rho, vphi, vz, Bphi соответственно
    np.savetxt("./out/data" + str(n) + ".dat", np.array(data).transpose(), fmt=('%.3e', '%.3e', '%.3e', '%.3e', '%.3e'))


print("Press Enter to start a simulation")
input()

n = 0 # номер шага по времени
SetIC()
# Сохранение НУ
SaveData(n)

n = n + 1
while t <= t_stop:
    contin = UpdateTimeStep()  
    # выход из цикла и сохранение результатов, если что-то пошло не так (TO DO)
    if(contin == False):
        print("Exiting")
        SaveData(n)
        break

    SetBC()
    
    Step()
    # вывод n, t, dt на экран и сохранение результатов каждые ... шагов
    if ((n % 10) == 0):
        print(  "step No ", n)
        print(  "t     = ", t, " s")
        print(  "dt    = ", dt, " s")
        print(  "v_max = ", v_max/1e5, " km/s")
        SaveData(n)
    UpdateIC()  
    t += dt
    n += 1

print("Total number of steps = ", n)
SaveData(n)