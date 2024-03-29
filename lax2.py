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
t_stop = 25
t = 0 #текущее время
c = 1.01 #Kurant number

#array for function values
un= np.zeros((N))
f= np.zeros((N))
un1 = np.zeros((N))
un_s = np.zeros((N))
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
        un[i] = U0(xs[i])

def SetBC():
    un1[0] = un[0] - 0.5*c*(un[1] - un[N-2])+0.5*c*c*(un[1]-2*un[0]+un[N-2])
    un1[N-1] = un[N-1] -0.5*c*(un[1] - un[N-2])+0.5*c*c*(un[1]-2*un[N-1]+un[N-2])



# расчет
    
def UpdateTimeStep():
    v_max = a_v
    dt = c*dx/v_max

def F(u):
    return a_v*u

def Step():
    for i in range ( 0, N-1 ):
        un_s[i] = 0.5*(un[i+1]+un[i])-(dt/dx)*0.5*(F(un[i+1])-F(un[i]))
    for i in range ( 1, N-1 ):
        un1[i] = un[i] - (dt/dx)*(F(un_s[i])-F(un_s[i-1]))
       
 
#обновление НУ
def UpdateIC():
    for i in range ( 0, N):
        un[i] = un1[i]


SetIC()
while t <= t_stop:
    SetBC()
    UpdateTimeStep()  
    Step()   
    UpdateIC() 
    t += dt
    print(t)

def SaveData():
    try:
        with open("data303.txt", "w") as f:
            f.write("#x u \n")
            for i in range(len(xs)):
                f.write(f"{xs[i]} {un1[i]} \n")
                #print(un1[i])
    except IOError:
        print("unable to open file for writing")
    #f.close()

SaveData()