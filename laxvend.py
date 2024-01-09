import numpy as np
import math 

import matplotlib.pyplot as plt 
import warnings
warnings.filterwarnings("ignore")

N=40
Nt = 10
h = 2*np.pi/N
k = 1/Nt
r = k/(h*h)
c = 0.9
time_steps=10
time=np.arange(0,(time_steps+.5)*k,k)
x=np.arange(0,2*np.pi+h/2,h)


X, Y = np.meshgrid(x, time)

fig = plt.figure()
plt.plot(X,Y,'ro');
plt.plot(x,0*x,'bo',label='Initial Condition');
plt.xlim((-h,2*np.pi+h))
plt.ylim((-k,max(time)+k))
plt.xlabel('x')
plt.ylabel('time (ms)')
plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
plt.title(r'Discrete Grid ',fontsize=24,y=1.08)
plt.show();

u=np.zeros((N+1,time_steps+1))
b=np.zeros(N-1)
# Initial Condition
for i in range (0,N+1):
    u[i,0]=1-np.cos(x[i])
    

fig = plt.figure(figsize=(8,4))
plt.plot(x,u[:,0],'o:',label='Initial Condition')
plt.xlim([-0.1,max(x)+h])
plt.title('Intitial Condition',fontsize=24)
plt.xlabel('x')
plt.ylabel('u')
plt.legend(loc='best')
plt.show()
ipos = np.zeros(N+1)
ineg = np.zeros(N+1)

for i in range(0,N+1):
   ipos[i] = int(i+1)
   ineg[i] = int(i-1)

ipos[N] = 0
ineg[0] = N


for j in range(0,time_steps):
    for i in range (0,N+1):
        u[i,j+1]=u[i,j]-c*(u[int(ipos[i]),j]-u[int(ineg[i]),j])+c*c/2*(u[int(ipos[i]),j]-2*u[i,j]+u[int(ineg[i]),j])
        
fig = plt.figure(figsize=(12,6))

plt.subplot(121)
for j in range (1,time_steps+1):
    plt.plot(x,u[:,j],'o:')
plt.xlabel('x')
plt.ylabel('u')

#plt.subplot(122)
#plt.imshow(u.transpose(), aspect='auto')
#plt.xticks(np.arange(len(x)), np.round(x,3),rotation=45)
#plt.yticks(np.arange(len(time)), np.round(time,2))
#plt.xlabel('x')
#plt.ylabel('time')
#clb=plt.colorbar()
#clb.set_label('w')
#plt.suptitle('Numerical Solution of the  Wave Equation'%(np.round(r,3)),fontsize=24,y=1.08)
#fig.tight_layout()
plt.show()