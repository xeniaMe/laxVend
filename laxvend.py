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
un= np.zeros((N))
un1 = np.zeros((N))
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
        un[i] = U0(xs[i])

def SetBC():
    un1[0] = un[0] - c*(un[0] - un[N-2])
    un1[N-1] = un[N-1] - c*(un[N-1] - un[N-2])


# расчет
    


def Step():
    for i in range ( 1, N-1 ):
        un1[i] = un[i] - c*(un[i+1]-un[i-1]) + (c*c/2)*(un[i+1] - 2*un[i]+un[i-1])
        print(un1[i])
 
        
        
        



#обновление НУ
def UpdateIC():
    for i in range ( 0, N):
        un[i] = un1[i]



while t <= t_stop:
    SetIC()
    SetBC()  
    Step()   
    UpdateIC() 
    t += dt

def SaveData():
    try:
        with open("C:\Рабочий стол\work1\data0_4_t.txt", "w") as f:
            f.write("#x u \n")
            for i in range(len(xs)):
                f.write(f"{xs[i]} {un1[i]} \n")
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