import numpy as np
import math 

import matplotlib.pyplot as plt 
import warnings
warnings.filterwarnings("ignore")


a = 0
b = 1
N = 80
v_max = 1
t_stop = 0.4
t = 0 #текущее время
c = 0.9 #Kurant number

#array for function values
un= np.zeros((N))
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
    un1[N-1] = un[N-1] -0.5*c*(un[1] - un[N-2])+0.5*c*c*(un[1]-un[N-1]+un[N-2])



# расчет
    


def Step():
    for i in range ( 0, N-1 ):
        un_s[i] = 0.5*(un[i+1]+un[i])-c*0.5*(un[i+1]-un[i])
        for i in range ( 1, N-1 ):
            un1[i] = un[i] - c*(un_s[i]-un_s[i-1])
       
 
#обновление НУ
def UpdateIC():
    for i in range ( 0, N):
        un[i] = un1[i]


SetIC()
while t <= t_stop:
    SetBC()  
    Step()   
    UpdateIC() 
    t += dt

def SaveData():
    try:
        with open("data300.txt", "w") as f:
            f.write("#x u \n")
            for i in range(len(xs)):
                f.write(f"{xs[i]} {un1[i]} \n")
                #print(un1[i])
    except IOError:
        print("unable to open file for writing")
    #f.close()

SaveData()