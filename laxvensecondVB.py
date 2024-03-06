import numpy as np
import math 

import matplotlib.pyplot as plt 
import warnings
warnings.filterwarnings("ignore")


a = 0
b = 1
N = 320
a_v= 0.7 #скорость переноса
v_max = a_v
t_stop = 0.4
t = 0 #текущее время
c = 1.01 #Kurant number
Bz =1
rho = 1

#array for function values
vn= np.zeros((N))
vn1 = np.zeros((N))
vn_s = np.zeros((N))
bn= np.zeros((N))
bn1 = np.zeros((N))
bn_s = np.zeros((N))
xs = np.linspace(a, b, N)


#построение расчетной сетки 
dx = (b-a)/ (N-1)
dt = c * dx/v_max
xs[0] = a

for i in range (1, N):
    xs[i] = xs[i-1] + dx


#определение начальных и граничных условий
def U0(x):
    if x <= 0.2:
        return 0.5
    elif x <= 0.4:
        return 1.0
    else:
        return 0.5


def SetIC():
    #global t, u, xs
    t = 0.0
    for i in range(N):
        vn[i] = U0(xs[i])
    for i in range(N):
        bn[i] = U0(xs[i])

def Fv(u):
    return -Bz/4*np.pi*rho*u

def Fb(u):
    return -Bz*u

def SetBC():
    vn1[0] = vn[0] - 0.5*c*(vn[1] - vn[N-2])+0.5*c*c*(vn[1]-2*vn[0]+vn[N-2])
    vn1[N-1] = vn[N-1] -0.5*c*(vn[1] - vn[N-2])+0.5*c*c*(vn[1]-2*vn[N-1]+vn[N-2])
    bn1[0] = bn[0] - 0.5*c*(bn[1] - bn[N-2])+0.5*c*c*(bn[1]-2*bn[0]+bn[N-2])
    bn1[N-1] = bn[N-1] -0.5*c*(bn[1] - bn[N-2])+0.5*c*c*(bn[1]-2*bn[N-1]+bn[N-2])



# расчет
    
def UpdateTimeStep():
    vA = Bz/np.sqrt(4*np.pi*rho)
    for j in range(N-1):
        v_max = max(abs(vn[j]) + vA)
    dt = c*dx/v_max



def Step1_v():
    for i in range ( 0, N-1 ):
        vn_s[i] = 0.5*(vn[i+1]+vn[i])-(dt/dx)*0.5*(Fv(vn[i+1])-Fv(vn[i]))
    for i in range ( 1, N-1 ):
        vn1[i] = vn[i] - (dt/dx)*(Fv(vn_s[i])-Fv(vn_s[i-1]))
       
 
def Step1_b():
    for i in range ( 0, N-1 ):
        bn_s[i] = 0.5*(bn[i+1]+bn[i])-(dt/dx)*0.5*(Fb(bn[i+1])-Fb(bn[i]))
    for i in range ( 1, N-1 ):
        bn1[i] = bn[i] - (dt/dx)*(Fb(bn_s[i])-Fb(bn_s[i-1]))

#обновление НУ
def UpdateIC():
    for i in range ( 0, N):
        vn[i] = vn1[i]
    for i in range ( 0, N):
        bn[i] = bn1[i]


SetIC()
while t <= t_stop:
    SetBC()
    UpdateTimeStep()  
    Step1_v()   
    Step1_b()   
    UpdateIC() 
    t += dt

def SaveData():
    try:
        with open("data305.txt", "w") as f:
            f.write("#x u \n")
            for i in range(len(xs)):
                f.write(f"{xs[i]} {vn1[i]} \n")
                #print(un1[i])
    except IOError:
        print("unable to open file for writing")
    #f.close()

SaveData()