"""
Date: April, 2019
Authors: Karina Chiñas Fuentes,
		 Daniel Hernández Mota

The code computes the energy and energy difference of N qubits
of a quantum computer based on a Carbon crystal with C13 particles.
The equations where derived using the Ising model.
"""

import numpy as np 
import matplotlib.pyplot as plt


# Number of qubits
n = int(input("Number of Qubits: "))

# Determine each frecuency of oscillation
w0 = 1/np.pi
w = []
for i in range(n):
    w.append(float(input("Frecuency for qubit " + str(n-i) + ": "))/w0)

# Determine the interaction
interaction=float((input("Interaction qubit-qubit: ")))/w0
print("\n")

J = []
for i in range(n-1):
    J.append(interaction/(10**i))
    
# Function for energy
def Energy(a, N):
    global w, J 
    a = np.array(a)
    w = np.array(w)
    one = (-np.ones(N))**a
    energy_NI = np.dot(one,w)
    a = list(a)
    energy_I = 0
    for j in range(N-1):
        for i in range(N-(j+1)):
            energy_I = energy_I + (J[j]*(-1)**(a[i]+a[i+(j+1)]))
        
    Energy_States =-0.5*( energy_NI + 0.5*energy_I )
    return Energy_States

# Function for generate all possibilities of bits (all states of qubits)
k = [ x for x in range(2**n) ]
k[0] = np.arange(1)
k[1] = np.arange(1,2)
m = 1
while m < n:
    for l in range(2**m):
        k[2**m + l] = np.concatenate(( np.arange(1), k[l] ))
        k[2**m + l][0] = 1 
        k[l] = np.concatenate(( np.arange(1), k[l] ))
    m += 1
a = [ list(x) for x in k ]

# Print the energies of the qubit. 
E=[]
for i in range(len(a)):
    E.append(Energy(a[i],n))
figure=[]   
for element in E:
    x=[.2,.8]
    y=[element,element]
    figure.append(plt.plot(x,y,"Black"))
    plt.xticks([])


# Energy plots 
title= "Energy States for "+str(n)+" Qubits"
plt.title(title)
plt.xlim(-.1,1)
fig = plt.gcf()
fig.set_size_inches(6,8)
plt.show()


# Interactions among the states of each qubits. 
vector = [ [] for i in range(2**n) ]

for b in range(n):
    c = b + 1
    for a in range(2**(c)):
        if a < 2**c/2:
            vector[a].append(int( (a + 1 + 2**(c-1) )%(2**c) ))
            if vector[a][b] == 0:
                vector[a][b] = 2**c
        elif a >= 2**(c-1):
            for i in range(len(vector[ int( a-2**(c-1) ) ])):
                vector[a].append( int( ( vector[int( a-2**(c-1) )][i] + 2**(c-1))%(2**c) ))
                if vector[a][i] == 0:
                    vector[a][i] = 2**c
    
vector = [sorted(x) for x in vector] 

states = input("Do you want to see the states interaction [Y/N]: ")

if (states == 'Y' or states == 'y' or states == 'Yes' or states == 'YES' or states == 'yes'):
    for i in range(len(vector)):
        print (i+1, ":" ,vector[i])


j = []
J = []

for i in range(2**n):
    for a in vector[i]:
        if (i+1) > a:
            diffE = E[i] - E[a-1]
            j.append( round(diffE,0) )
            J.append([i+1,a]) 


# Allow us to know which energy differences are equal among qubits 
s = 0
while s < len(j): 
    for l in range(s,len(j)):
        
        if s == l:
            pass

        elif j[l] == j[s]: 
            print( "WARNING: The energy difference " + str(J[s]) + " is equal to " + str(J[l]) )

    s += 1

print("\n")
for i in range(len(J)):
    print(j[i],J[i])

# If you want to see all transitions, uncomment the following lines:
#        if (i+1) < a:
#            print ("E [",a,"]-E [",i+1,"]")