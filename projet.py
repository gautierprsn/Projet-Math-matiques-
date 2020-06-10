import numpy as np
import scipy.integrate as integr
import matplotlib.pyplot as plt
import scipy.stats as stats
import random
from mpl_toolkits.mplot3d import Axes3D
from math import *


def binom(n,i):
    if n<i: return False
    I=1
    N=1
    S=1
    for k in range(1,i+1):
        I*=k
    for j in range(1,n+1):
        N=N*j
    for t in range(1,n-i+1):
        S*=t
    return N/(S*I)
    
def f(x):
    return max(100-x,0)

## Question 3
def pricer1(rN,N,bN,hN,s,f):
    somme=0
    qN=(rN-bN)/(hN-bN)
    for k in range(0,30):
        somme+=binom(N,k)*f(s*(1+hN)**k*(1+bN)**(N-k))*qN**k*(1-qN)**(N-k)
    return  1/(1+rN)**N*somme
    
    
## Question 5
def pricer2(f,N,rN,bN,hN,s):
    Vn=np.eye(N+1,N+1)
    qN=np.array((rN-bN)/(hN-bN))
    for k in range (0,N+1):
        Vn[k,N]=f(s*(1+bN)**(k)*(1+hN)**(N-k))
    for i in range(N,-1,-1):
        for j in range(N-1,-1,-1):
            if (i>j):
                Vn[i,j]=0
            else:
                Vn[i,j]=np.array(1/(1+rN))*[qN*Vn[i,j+1] + (np.array(1)-qN)*Vn[i+1,j+1]]
    return Vn[0,0]

## Question 7
def compare (f,rN,bN,hN,s):
    N=random.randint(5,15)
    print (N);
    return pricer2(f,N,rN,bN,hN,s)- pricer1(rN,N,bN,hN,s,f)
    
## Question 9 
def couverture(N,s,rN,hN,bN,f):
    Sh = pricer2(f,N,rN,bN,hN,(1+hN)*s);
    Sb = pricer2(f,N,rN,bN,hN,(1+bN)*s);
    alpha =  (Sh-Sb)/(s*(hN-bN));
    beta = (Sh*(1+bN)-Sb*(1+hN))/((1+rN)**N *(bN-hN));
    print(alpha)
    print(beta)

    
    
## Question 12
def pricer3(n,r,s,T,sigma,f):
    somme=0
    for i in range(1,n):
        somme=somme+exp(-r*T)*f(s*exp((r-(sigma**2)/2)*T+sigma*sqrt(T)*stats.norm.rvs(0,1)))
    return (1/n)*somme
    

    
## Question 15
def put(s,r,sigma,T,K):
    d=(1/sigma*sqrt(T))*(log(s/K)+((r+sigma**2)/2)*T)
    return -s*stats.norm.cdf(-d, loc = 0, scale = 1)+K*exp(-r*T)*stats.norm.cdf(-d+sigma*sqrt(T), loc = 0, scale = 1)

## Question 13 et 17

for k in range (1,10):
    plt.plot((10**5)*k, pricer3((10**5)*k,0.01,100,1,0.1,f), marker='*', linestyle='solid', color='r')
plt.axhline(y=put(100,0.01,0.1,1,100))
plt.xlabel('N')
plt.ylabel('Pricer3 ')

plt.title('Comparaison du pricer3 avec le put en fonction de N ')
plt.show()

## Question 18

t=np.array([1/12,1/6,1/4,1/3,1/2,1])
s=np.array([20,40,60,80,100,120,140,160,180,200])
fig = plt.figure()
 
ax = fig.gca(projection='3d')
Put=np.array([0 for k in range (1,61)])
for i in range(0,10) :
    for j in range (0,6):
        ax.scatter(t[j],s[i],put(s[i],0.01,0.1,t[j],100),cmap='hot')
ax.set_xlabel('T')
ax.set_ylabel('s')
ax.set_zlabel('put')
 
plt.show()


## Question 19

# def fT(t):
#     s=100
#     sigma=0.3
#     T=1
#     r=0.02
#     St=s*exp((r-(sigma**2)/2)*T+sigma*sqrt(T)*stats.norm.rvs(loc=0,scale=1,size=1))
#     return max(100-St,0)
# 
# for k in range (1,100):
#     N=10*k
#     s=100
#     sigma=0.3
#     r=0.02
#     T=1
#     rN=r*T/N
#     hN=(1+rN)*exp(sigma*sqrt(T/N))-1
#     bN=(1+rN)*exp(-sigma*sqrt(T/N))-1
#     plt.plot(N, pricer2(fT,10*k,rN,bN,hN,s), marker='*',color='r')
# plt.axhline(y=put(100,r,sigma,1,100), marker='*')
# plt.xlabel('N')
# plt.ylabel('pricer2')
# 
# plt.title('pricer2 en fonction de N ')
# plt.show()



 

