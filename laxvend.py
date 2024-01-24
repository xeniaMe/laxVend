import numpy as np
import math 

import matplotlib.pyplot as plt 
import warnings
warnings.filterwarnings("ignore")


a = 0
b = 1
N = 40
v = 1
t_stop = 0.4
t = 0 #текущее время
c = 0.9 #Kurant number

#array for function values
u= np.zeros((N))
u1 = np.zeros((N))
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
        u[i] = U0(xs[i])

def SetBC():
    u1[0] = 5 #u[0] - c*(u[0] - u[N-2])
    u1[N-1] = 10 #u[N-1] - c*(u[N-1] - u[N-2])


# расчет
def Step():
    for i in range ( 1, N-1 ):
        u1[i] = u[i] - c*(u[i+1]-u[i-1]) + c*c/2*(u[i+1] - 2*u[i]+u[i-1])
    #print(u1)
        

#обновление НУ
def UpdateIC():
    for i in range ( 0, N):
        u[i] = u1[i]



while t <= t_stop:
    SetBC()  
    Step()   
    UpdateIC()  
    t += dt

def SaveData():
    try:
        with open("data.txt", "w") as f:
            f.write("#x u\n")
            for i in range(len(xs)):
                f.write(f"{xs[i]} {u1[i]}\n")
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