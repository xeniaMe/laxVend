import numpy as np
import math 

import matplotlib.pyplot as plt 
import warnings
warnings.filterwarnings("ignore")


L = 1e6 # длина расчетной области 
a = 0   # координата левой границы расчетной области
b = L   # координата правой границы расчетной области
N = 320
a_v= 0.7 #скорость переноса
v_max = a_v
t = 0 #текущее время
c = 0.95 #Kurant number
Bz = 0.5 #Г, МП
rho = 6*1.672e-24 #плотность протонов

vA = Bz/np.sqrt(4*np.pi*rho)
t_stop = L/vA
v_0 = 0.05*vA
lambd = 0.1*L
k = 2*np.pi/lambd

print("Model parameters:")
print(" B0     = ", Bz, " G")
print(" rho0   = ", rho, " g/cm^3")
print(" L      = ", L / 1e5, " km")
print(" vA     = ", vA / 1e5, " km/s")
print(" t_stop = ", t_stop, " s")
print(" lambda = ", lambd / 1e5, " km")

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
xs[0] = a

for i in range (1, N):
    xs[i] = xs[i-1] + dx


#определение начальных и граничных условий
    
def b0(x):
    return 0.0  # Начальное Bx пусть будет 0

def v0(x):
    return  v_0 * np.sin(x * k)
    


def SetIC():
    #global t, u, xs
    t = 0.0
    for i in range(N):
        vn[i] = v0(xs[i])
        bn[i] = b0(xs[i])

def Fv(b):
    return -Bz/(4*np.pi*rho)*b

def Fb(u):
    return -Bz*u

def SetBC():
    vn1[0] = vn[0] - 0.5*c*(vn[1] - vn[N-2])+0.5*c*c*(vn[1]-2*vn[0]+vn[N-2])
    vn1[N-1] = vn[N-1] -0.5*c*(vn[1] - vn[N-2])+0.5*c*c*(vn[1]-2*vn[N-1]+vn[N-2])
    bn1[0] = bn[0] - 0.5*c*(bn[1] - bn[N-2])+0.5*c*c*(bn[1]-2*bn[0]+bn[N-2])
    bn1[N-1] = bn[N-1] -0.5*c*(bn[1] - bn[N-2])+0.5*c*c*(bn[1]-2*bn[N-1]+bn[N-2])



# расчет
    
def UpdateTimeStep():
    for j in range(N-1):
        v[j] = abs(vn[j]) + vA
    v_max = max(v)
    dt = c*dx/v_max



def Step():
    for i in range ( 0, N-1 ):
        vn_s[i] = 0.5*(vn[i+1]+vn[i])-(dt/dx)*0.5*(Fv(bn[i+1])-Fv(bn[i]))
        bn_s[i] = 0.5*(bn[i+1]+bn[i])-(dt/dx)*0.5*(Fb(vn[i+1])-Fb(vn[i]))    
    for i in range ( 1, N-1 ):
        vn1[i] = vn[i] - (dt/dx)*(Fv(bn_s[i])-Fv(bn_s[i-1]))  
        bn1[i] = bn[i] - (dt/dx)*(Fb(vn_s[i])-Fb(vn_s[i-1]))
    # print(vn)

#обновление НУ
def UpdateIC():
    for i in range ( 0, N):
        vn[i] = vn1[i]
        bn[i] = bn1[i]


SetIC()
while t <= t_stop:
    #print(  "t = ", t, " s")
    #print(  "  dt = ", dt, " s")
    SetBC()
    UpdateTimeStep()  
    Step()     
    UpdateIC()  
    t += dt

def SaveData():
    try:
        with open("data305.txt", "w") as f:
            f.write("#x u \n")
            for i in range(len(xs)):
                f.write(f"{xs[i]} {vn1[i]} \n")
    except IOError:
        print("unable to open file for writing")
    #f.close()

SaveData()