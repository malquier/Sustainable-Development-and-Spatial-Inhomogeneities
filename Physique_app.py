##Bibliothèque :
    
import numpy as np
import matplotlib.pyplot as plt
import random

## Initialisation :
mu = 0.3
K1 = 100
rho = 0.2
phi = 10
P0 = 10
deltA = 0.03
deltK = 0.15
deltR = 0.3 
f = 0.5
noise = 0.1
p = 10
f1 = 0.2
pas = 1
d12 = 1 
t0, tf = 0, 10000
N = int((tf-t0)/pas)
time = np.linspace(t0,tf,N)
K0 = phi*mu/((1+deltR)*deltK*deltA) - K1 
R0 = phi/(1+deltR)
A0 = abs(mu/deltA * (1 - K1*(1+deltR)*deltK/phi))
[K,A,R,C] = np.zeros(N),np.zeros(N),np.zeros(N),np.zeros(N)


## The AB model

""" We first write functions useful for the resolution of the differential equations"""

def C_opt(K,A,p):
    return K*(np.sqrt(A/p)-1)

def prod(A,K,C):
    return (A*K*C/(K+C))
    
def profit(A,K,C,p):
    return(prod(A,K,C)-p*C)

def attract(f,d12):
    return(f**d12/(f**d12+f**d12))
 
""" We initialize the differents variables used as capital, consumption..."""   

K[0] = np.random.poisson(K0)
A[0] = np.random.poisson(A0)
R[0] = np.random.poisson(R0)
C0 = min(C_opt(K0,A0,p) , R[0])
C[0] = C0
 
""" We can start the resolution of the differential equation """

for i in range(N-1):
    # Noise :
    mu=np.random.normal(mu,noise*mu)
    rho=np.random.normal(rho,noise*rho)
    
    # Méthode d'Euler
    K[i+1] = K[i] + pas * rho * (prod(A[i],K[i],C[i]) - p*C[i] - deltK*K[i])
    A[i+1] = A[i] + pas * (mu*K[i]/(K[i]+K1) - deltA*A[i])
    R[i+1] = R[i] + pas * (phi - C[i] - deltR*R[i])
    C[i+1] = min(C_opt(K[i+1],A[i+1],p) , R[i+1])
    print(K[i+1],C[i],R[i])

plt.plot(time,(K))
plt.plot(time,(A))
plt.grid(True)
plt.xlabel("Temps")
plt.ylabel("Capital et coefficient proportionnel à la technologie")
plt.title("Evolution du capital et de la technicité en fonction du temps")
plt.show()
    

## The spatial model :
    
"""We use the same model AK but with two local economies : they are exportations, importations"""
[K3,A1,R1,C1] = np.zeros(N),np.zeros(N),np.zeros(N),np.zeros(N)
[K2,A2,R2,C2] = np.zeros(N),np.zeros(N),np.zeros(N),np.zeros(N)
K3[0] = np.random.poisson(K0)
A1[0] = np.random.poisson(A0)
R1[0] = np.random.poisson(R0)
K2[0] = np.random.poisson(K0)
A2[0] = np.random.poisson(A0)
R2[0] = np.random.poisson(R0)
C0 = min(C_opt(K[0],A[0],p) , R[0])
C1[0] = C0
C12, C21 = 0,0

for i in range(N-1):
    K3[i+1] = K3[i] + pas * rho * (prod(A1[i],K3[i],C1[i]-C12+C21) - p*(C1[i]-C12+C21) - deltK*K3[i])
    A1[i+1] = A1[i] + pas * (mu*K3[i]/(K3[i]+K1) - deltA*A1[i])
    R1[i+1] = R1[i] + pas * (phi - C1[i] - deltR*R1[i])
    K2[i+1] = K2[i] + pas * rho * (prod(A2[i],K2[i],C2[i]-C21+C12) - p*(C2[i]-C21+C12) - deltK*K2[i])
    A2[i+1] = A2[i] + pas * (mu*K2[i]/(K2[i]+K1) - deltA*A2[i])
    R2[i+1] = R2[i] + pas * (phi - C2[i] - deltR*R2[i])
    C1[i+1] = min(C_opt(K3[i+1],A1[i+1],p) , R1[i+1])
    C2[i+1] = min(C_opt(K2[i+1],A2[i+1],p) , R2[i+1])
    C12 = attract(f,d12)*C1[i]
    C21 = attract(f,d12)*C2[i]
    
    
# plt.plot(time,A1)
# plt.plot(time,A2)
# plt.grid(True)
# plt.ylabel("Capital des villes 1 et 2")
# plt.xlabel("Temps")
# plt.title("Evolution du capital de deux villes en fonction du temps")
# plt.show()