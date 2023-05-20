

import numpy as np
import math
import time
import Solution_test

MAX = 100

Ans=[]
Ans.clear()







def crout(A,b,p):
        print("\n Crout Method Steps : ")
        Ans.clear()
        unique, infinite = Solution_test.IsUnique(A, b)
        if (not unique and not infinite):
         x = []
         flag = "System is INCONSISTENT"
         x.append(flag)
         return x, None
        if (infinite):
         x = []
         flag = "System has Infinite Solution"
         x.append(flag)
         return x, None
        a=np.array(A)
        L = np.zeros((len(a), len(a)))
        U = np.zeros((len(a), len(a)))
        t1=time.time()
        for k in range(0, len(a)):
            U[k, k] = 1 

            for j in range(k, len(a)):
                sum0 = sum([L[j, s] * U[s, k] for s in range(0, j)]) #range from index 0
                L[j, k] = a[j, k] - sum0 #reversed index
            for j in range(k+1, len(a)):
                sum1 = sum([L[k, s] * U[s, j] for s in range(0, k)]) #range from index 0
                U[k, j] = (a[k, j] - sum1) / L[k, k]
        print("Upper Triangular Matrix")
        print(U)
        print("\n")
        print("Lower Triangular Matrix")
        print(L)
        x=solve(L,U,b,p)
        t2 = time.time()
        x.append("correct")
        return x,t2-t1

    


def solve(L,U,b,p):
    start_time = time.time()
    L = np.array(L)
    U = np.array(U)
    b = np.array(b).transpose()
    print("L=")
    print(L)
    print("\n")
    print("U=")
    print(U)
    print("\n")
    print("b=",end="")
    print(b)
    print("\n")

    
    y = [0 for _ in range(len(L))]
    i = 1
    print("Y=",end="")
    print(y)
    print("\n")
    y[0] = round((b[0]/L[0][0]),p)
    print("Y=",end="")
    print(y)
    print("\n")
    
    for i in range(1,len(y),1):
        sum0 = 0
        j=0
        for j in range(0,i,1):
            sum0 += round((L[i][j]*y[j]),p)
            y[i] = round(((b[i] - sum0)/L[i][i]),p)
            y = np.round_(y,decimals=p)
            print("Y= ",end="")
            print(y)

        
        
        
        


    print("Y= ",end="")    
    print(y)
    print("\n")
    x = [0 for _ in range(len(y))]
    
    x[len(y)-1] = round((y[len(y)-1]/U[len(y)-1][len(y)-1]),p)
    print("X= ",end="")
    print(x)
    for j in range(len(y)-2,-1,-1):
        sum1 = 0
        for q in range(j+1,len(y),1):
            sum1 += round((U[j][q]*x[q]),p)
            x[j] = round((y[j] - sum1)/U[j][j],p)
            x = np.round_(x,decimals=p)
            print("X= ",end="")
            print(x)
    print("X= ",end="")
    print(x)
    for i in range(len(x)):
        Ans.append(x[i])
    return Ans


    

