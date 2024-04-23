import numpy as np
import math 

import matplotlib.pyplot as plt 
import warnings
warnings.filterwarnings("ignore")


L = 1e6 # длина расчетной области , см
a = 0   # координата левой границы расчетной области
b = L   # координата правой границы расчетной области
N = 320

t = 0 #текущее время
c = 0.5 #Kurant number
Bz = -1 #Гс, минус - для того, чтобы волна распространялась "вправо"
n = 1e8 # концентрация частиц 
rho = n*1.672e-24 # плотность протонов
print(rho)

vA = np.abs(Bz)/np.sqrt(4*np.pi*rho)
v_0 = 0.05*vA
B_0 = 0.05*np.abs(Bz)
t_A = L/vA
t_stop = 4*t_A 
lambd = 0.1*L
k = 2*np.pi/lambd
# начальная максимальная скорость на сетке
v_max = abs(vA)

print("Model parameters:")
print(" B0     = ", Bz, " G")
print(" rho0   = ", rho, " g/cm^3")
print(" v0     = ", v_0 / 1e5, " km/s")
print(" L      = ", L / 1e5, " km")
print(" vA     = ", vA / 1e5, " km/s")
print(" t_stop = ", t_stop, " s")
print(" lambda = ", lambd / 1e5, " km")
print(" v_max  = ", v_max / 1e5, " km/s")

#array for function values
vn= np.zeros((N))
vn1 = np.zeros((N))
vn_s = np.zeros((N))
v= np.zeros((N))
bn= np.zeros((N))
bn1 = np.zeros((N))
bn_s = np.zeros((N))
xs = np.linspace(a, b, N)


#построение расчетной сетки 
dx = (b-a)/ (N-1)
dt = c*dx/v_max
print(" dx = ", dx / 1e5, " km")
print(" dt = ", dt, " s")

xs[0] = a
for i in range (1, N):
    xs[i] = xs[i-1] + dx
    
#определение начальных и граничных условий    
def b0(x):
    return B_0 * np.sin(x * k) # Начальное Bx пусть будет 0

def v0(x):
    return  v_0 * np.sin(x * k)
    
x0 = 0.0#0.5*L - lambd
def SetIC():
    global t, vn, bn, xs
    t = 0.0
    for i in range(N):
        # задается одно колебание в области [x0, x0+lambd]
        if ((xs[i] > x0) and (xs[i] < lambd+x0)):
            vn[i] = v0(xs[i])
            bn[i] = b0(xs[i])
        else:
            vn[i] = 0.0
            bn[i] = 0.0
            

def Fv(b):
    return (-Bz/(4*np.pi*rho))*b

def Fb(u):
    return -Bz*u

def SetBC():
    global t, vn_s, bn_s, vn1, bn1, vn, bn, dt, dx
    # периодические ГУ
    # левая граница
    vn_s[N-2] = 0.5*(vn[N-1]+vn[N-2])-(dt/dx)*0.5*(Fv(bn[N-1])-Fv(bn[N-2]))
    bn_s[N-2] = 0.5*(bn[N-1]+bn[N-1])-(dt/dx)*0.5*(Fb(vn[N-1])-Fb(vn[N-2]))
    vn_s[0] = 0.5*(vn[1]+vn[0])-(dt/dx)*0.5*(Fv(bn[1])-Fv(bn[0]))
    bn_s[0] = 0.5*(bn[1]+bn[0])-(dt/dx)*0.5*(Fb(vn[1])-Fb(vn[0]))    
    vn1[0] = vn[0] - (dt/dx)*(Fv(bn_s[0])-Fv(bn_s[N-2]))  
    bn1[0] = bn[0] - (dt/dx)*(Fb(vn_s[0])-Fb(vn_s[N-2]))
    # правая граница
    vn_s[N-1] = 0.5*(vn[1]+vn[N-1])-(dt/dx)*0.5*(Fv(bn[1])-Fv(bn[N-1]))
    bn_s[N-1] = 0.5*(bn[1]+bn[N-1])-(dt/dx)*0.5*(Fb(vn[1])-Fb(vn[N-1]))    
    vn1[N-1] = vn[N-1] - (dt/dx)*(Fv(bn_s[N-1])-Fv(bn_s[N-2]))  
    bn1[N-1] = bn[N-1] - (dt/dx)*(Fb(vn_s[N-1])-Fb(vn_s[N-2]))

# расчет
def UpdateTimeStep():
    global dt, v_max, v, vn
    for j in range(N-1):
        v[j] = abs(vn[j]) + vA
    v_max = max(v)
    dt = c*dx/v_max
    # если макс. скорость становится слишком большой - код нужно остановить
    if v_max > 5*vA:
        return False
    else:
        return True

def Step():
    global dt, vn_s, bn_s, vn1, bn1, vn, bn, dt, dx
    for i in range ( 0, N-1 ):
        vn_s[i] = 0.5*(vn[i+1]+vn[i])-(dt/dx)*0.5*(Fv(bn[i+1])-Fv(bn[i]))
        bn_s[i] = 0.5*(bn[i+1]+bn[i])-(dt/dx)*0.5*(Fb(vn[i+1])-Fb(vn[i]))    
    for i in range ( 1, N-1 ):
        vn1[i] = vn[i] - (dt/dx)*(Fv(bn_s[i])-Fv(bn_s[i-1]))  
        bn1[i] = bn[i] - (dt/dx)*(Fb(vn_s[i])-Fb(vn_s[i-1]))
    # print(vn)

#обновление НУ
def UpdateIC():
    global t, vn, bn, vn1, bn1
    for i in range ( 0, N):
        vn[i] = vn1[i]
        bn[i] = bn1[i]

def SaveData(n):
    data = []
    data.append(xs)

    if (n == 0): # начальное условие
        data.append(vn)
        data.append(bn)
    else:
        data.append(vn1)
        data.append(bn1)
        
    # данные сохраняются в папку ./out/
    # столбцы 1, 2, 3 - массивы xs, vn1, bn1 соответственно
    np.savetxt("./out/data_v" + str(n) + ".dat", np.array(data).transpose(), fmt=('%.3e', '%.3e', '%.3e'))


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
    if ((n % 20) == 0):
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