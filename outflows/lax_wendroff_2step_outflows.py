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
N = 320
t = 0 #текущее время
c = 0.5 #Kurant number

L = 5   # длина расчетной области , H
a = 0   # координата левой границы расчетной области
b = L   # координата правой границы расчетной области

# ----- параметры физ. модели -------------------------------------------------
r_cm = 1 * au # радиальное расстояние от звезды
rho_ISM = np.exp(-0.5) # безразмерная плотность МЗС, в ед. rho0
M_k = 1.0
beta = 10.0 # плазменный параметр
Bz = 1.0

# ----- масштабы безразмерных переменных --------------------------------------
rho0 = 1e-13
T0 = 500
c_T = np.sqrt(Rg * T0 / mu)
v0 = c_T
# L0 = rho0 * r_cm * v0
H = 0.05 * r_cm
t0 = H / v0
B0 = 1
g0 = v0**2/H

# ----- коэффициенты уравнений ------------------------------------------------
def v_k(r_au):
    return np.sqrt(G * Mstar / (r_au * au))

vk0 = v_k(r_cm / au)

def v_phi0(r_cm, z_cm):
    return v_k(r_cm / au) * np.power(1.0 + (z_cm / r_cm)**2, -0.75)

def v_phi(L, rho):
    return L / rho

def rho_hs(z):
    return np.exp(-0.5 * z**2)

def p(rho):
    return rho

def g_z(r_cm, z_cm, M_g):
    return G*M_g*z_cm/np.power(r_cm**2 + z_cm**2, 1.5)

vA = B0 / np.sqrt(4.0 * np.pi * rho0)
t_stop = 0.01
# начальная максимальная скорость на сетке
v_max = abs(vA)

print("Model scales:")
print(" B0     = ", B0, " G")
print(" rho0   = ", rho0, " g/cm^3")
# print(" L0      = ", L0, " g/cm/s")
print(" t0 = ", t0/year, " yr")
print(" v0 = ", v0/1e5, " km/s")

print("Model charcteristics:")
print(" vA     = ", vA / 1e5, " km/s")
print(" v_k    = ", vk0 / 1e5, " km/s")
print(" M_k    = ", vk0 / c_T)
print(" t_stop = ", t_stop, " [t0]")
print(" v_max  = ", v_max / 1e5, " km/s")

# ----- основные переменные ---------------------------------------------------

# консервативные переменные на шаге t^n
u0_n = np.zeros((N)) # = rho
u1_n = np.zeros((N)) # = rho*v_z
u2_n = np.zeros((N)) # = rho*v_phi
u3_n = np.zeros((N)) # = bphi
# консервативные переменные на шаге t^(n+1/2) (промежуточном)
u0_s = np.zeros((N)) # = rho
u1_s = np.zeros((N)) # = rho*v_z
u2_s = np.zeros((N)) # = rho*v_phi
u3_s = np.zeros((N)) # = bphi
# консервативные переменные на шаге t^(n+1)
u0_n1 = np.zeros((N)) # = rho
u1_n1 = np.zeros((N)) # = rho*v_z
u2_n1 = np.zeros((N)) # = rho*v_phi
u3_n1 = np.zeros((N)) # bphi

# списки векторов консервативных переменных
u_n  = [u0_n, u1_n, u2_n, u3_n]
u_s  = [u0_s, u1_s, u2_s, u3_s]
u_n1 = [u0_n1, u1_n1, u2_n1, u3_n1]

# примитивные переменные на шаге t^n
rho_n  = np.zeros((N))  # = rho
vz_n    = np.zeros((N)) # = v_z
vphi_n = np.zeros((N))  # v_phi
B_n    = np.zeros((N))  # bphi

# rho_s = np.zeros((N)) # = rho
# v_s = np.zeros((N)) # =v_z
# L_s = np.zeros((N)) # L = rho*r*v_phi
# vphi_s = np.zeros((N)) # v_phi
# B_s = np.zeros((N)) # bphi

# rho_n1 = np.zeros((N)) # = rho
# v_n1 = np.zeros((N)) # = v_z
# L_n1 = np.zeros((N)) # L = rho*r*v_phi
# vphi_n1 = np.zeros((N)) # v_phi
# B_n1 = np.zeros((N)) # bphi

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
def rho_IC(z):
    if z < 1.0:
        return np.exp(-0.5 * z**2)
    else:
        return np.exp(-0.5)

def v_IC(z):
    return 0.0

def vphi_IC(z):
    if z < 1.0:
        return v_phi0(r_cm, z * H) / v0
    else:
        return 0.0
    
# def L_IC(z):
#     if z < 1.0:
#         return rho_IC(z) * r_cm * v_phi0(r_cm, z * H) / L0
#     else:
#         return 0.0
    
def B_IC(z):
    return 0.0

def p2c(var_n, cv, pv):
    """
    Перевод вектора примитивных переменных в векор консервативных переменных для переменной с номером var_n = {0, 1, 2, 3}

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
    # pv = {rho, vz, vphi, Bphi}
    
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
     Перевод вектора консервативных переменных в векор примитивных переменных для переменной с номером var_n = {0, 1, 2, 3}

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
    for i in range(N):
        pv[0][i] = rho_IC(zs[i])  # rho_n
        pv[1][i] = vphi_IC(zs[i]) # vphi_n
        pv[2][i] = v_IC(zs[i])    # vz_n
        pv[3][i] = B_IC(zs[i])    # B_n
        
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
        return -g_z(r_cm, zs[i]*H, Mstar)/g0 * u[0][i]
    elif var_n == 3:
        return 0.0
    else:
        None
        
def SetBC():
    global t, u_n1

    # левая граница
    rho_L  = 1.0 # 
    vphi_L = vk0 / v0 #
    vz_L    = 0.0 #
    B_L    = 0.0 #
    
    u_n1[0][0] = rho_L
    u_n1[1][0] = rho_L * vphi_L
    u_n1[2][0] = rho_L * vz_L
    u_n1[3][0] = B_L
    
    # правая граница
    rho_R  = rho_ISM
    vphi_R = 0.0 #
    vz_R   = 0.0 #
    B_R    = 0.0 #
    
    u_n1[0][N-1] = rho_R
    u_n1[1][N-1] = rho_R * vphi_R
    u_n1[2][N-1] = rho_R * vz_R
    u_n1[3][N-1] = B_R
    

# расчет
def UpdateTimeStep():
    global dt, v_max, v, v_n

    for i in range(N):
        # квадрат безразмерной альв. скорости
        vAvA = (2.0/beta) * (Bz**2 + B_n[i]**2) / rho_n[i]
        # полная скорость, 1.0 - это б.р. скорость звука
        vv[i] = abs(vz_n[i]) + abs(vphi_n[i]) + np.sqrt(1.0 + vAvA)
        
    v_max = max(vv)
    dt = c*dz/v_max
    # # если макс. скорость становится слишком большой - код нужно остановить
    # if v_max > 5*vA:
    #     return False
    # else:
    #     return True

def Step():
    global dt, u_n, u_s, u_n1, dz
    
    for var_n in range(4):
        for i in range(0, N-1):
            u_s[var_n][i] =  0.5 * (u_n[var_n][i+1] + u_n[var_n][i]) - 0.5*(dt/dz)*(F(u_n, i+1, var_n) - F(u_n, i, var_n)) + Source(u_n, i, var_n)
            
    for var_n in range(4):
        for i in range (1, N-1):
            u_n1[var_n][i] = u_n[var_n][i] - (dt/dz) * (F(u_s, i, var_n) - F(u_s, i-1, var_n)) + 0.5 * (Source(u_s, i, var_n) + Source(u_n, i, var_n))
            
    for var_n in range(4):
        c2p(var_n, u_n1, pv)

#обновление НУ
def UpdateIC():
    global t, u_n, u_n1
    for var_n in range(4):
        u_n[var_n] = u_n1[var_n]

def SaveData(n):
    data = []
    data.append(zs)
    
    if (n == 0): # начальное условие
        for var_n in range(4):
            c2p(var_n, u_n, pv)
            data.append(pv[var_n])
    else:
        for var_n in range(4):
            c2p(var_n, u_n1, pv)
            data.append(pv[var_n])
        
    # данные сохраняются в папку ./out/
    # столбцы 1, 2, 3 - массивы xs, vn1, bn1 соответственно
    np.savetxt("data" + str(n) + ".dat", np.array(data).transpose(), fmt=('%.3e', '%.3e', '%.3e', '%.3e', '%.3e'))


print("Press Enter to star a simulation")
input()

n = 0 # номер шага по времени
SetIC()
# Сохранение НУ
SaveData(n)

n = n + 1
while t <= t_stop:
    contin = UpdateTimeStep()  
    # выход из цикла и сохранение результатов, если скорость стала слишком большой
    if(contin == False):
        print("Alfven speed became too large, exiting")
        SaveData(n)
        break

    SetBC()
    
    Step()
    # вывод n, t, dt на экран и сохранение результатов каждые ... шагов
    if ((n % 1) == 0):
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