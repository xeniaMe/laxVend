import numpy as np
import math 

import matplotlib.pyplot as plt 
import warnings
warnings.filterwarnings("ignore")


a = 0
b = 1
N = 80
v = 1
t_stop = 0.4
t = 0 #текущее время
c = 0.9 #Kurant number

#array for function values
Vn= np.zeros((N))
Bn= np.zeros((N))
#un0_5= np.zeros((N))
#bn= np.zeros((N))
V_new = np.zeros((N+1))
B_new = np.zeros((N + 1))
xs = np.linspace(a, b, N)


#построение расчетной сетки 
dx = (b-a)/ (N-1)
dt = c * dx/v
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
        Vn[i] = U0(xs[i])
        Bn[i] = U0(xs[i])


def SetBC():
    V_new[0] = Vn[0] - 0.5*c*(Vn[1] - Vn[N-2])+0.5*c*c*(Vn[1]-2*Vn[0]+Vn[N-2])
    V_new[N-1] = Vn[N-1] -0.5*c*(Vn[1] - Vn[N-2])+0.5*c*c*(Vn[1]-Vn[N-1]+Vn[N-2])
    B_new[0] = Bn[0] - 0.5*c*(Bn[1] - Bn[N-2])+0.5*c*c*(Bn[1]-2*Bn[0]+Bn[N-2])
    B_new[N-1] = Bn[N-1] -0.5*c*(Bn[1] - Bn[N-2])+0.5*c*c*(Bn[1]-Bn[N-1]+Bn[N-2])



# расчет
    


def Step1():
    #for n in range(N):
    for i in range(1, N-1):
        V_new[i] = 0.5 * (Vn[i + 1] + Vn[i - 1]) - 0.5 * a * dt / dx * (Bn[i + 1] - Bn[i - 1])
        B_new[i] = 0.5 * (Bn[i + 1] + Bn[i - 1]) - 0.5 * b * dt / dx * (Vn[i + 1] - Vn[i - 1])
    Vn[:] = V_new
    Bn[:] = B_new
#обновление НУ
def UpdateIC():
    for i in range ( 0, N):
        Vn[i] = V_new[i]
        Bn[i] = B_new[i]


SetIC()
while t <= t_stop:
    SetBC()  
    Step1()   
    UpdateIC() 
    t += dt

def SaveData():
    try:
        with open("data100.txt", "w") as f:
            f.write("#x u \n")
            for i in range(len(xs)):
                f.write(f"{xs[i]} {V_new[i]} \n")
                #print(un1[i])
    except IOError:
        print("unable to open file for writing")
    #f.close()

SaveData()






#plt.plot(x,0*x,'bo',label='Initial Condition');
# plt.xlim((-h,2*np.pi+h))
# plt.ylim((-k,max(time)+k))
# plt.xlabel('x')
# plt.ylabel('time (ms)')
# plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
# plt.title(r'Discrete Grid ',fontsize=24,y=1.08)
# plt.show();
# fig = plt.figure(figsize=(8,4))
# plt.plot(x,u[:,0],'o:',label='Initial Condition')
# plt.xlim([-0.1,max(x)+h])
# plt.title('Intitial Condition',fontsize=24)
# plt.xlabel('x')
# plt.ylabel('u')
# plt.legend(loc='best')
# plt.show()       
# fig = plt.figure(figsize=(12,6))
# plt.subplot(121)
# for j in range (1,time_steps+1):
#     plt.plot(x,u[:,j],'o:')
# plt.xlabel('x')
# plt.ylabel('u')
#