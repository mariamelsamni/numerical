

import numpy as np
import math
import Solution_test
import time

MAX = 100

def Hermitian_Positive_Semidefinite(matrix,n):
    flag = True
    flag1 = True
    
    eigen_values, eigen_vectors = np.linalg.eig(matrix)
    for j in range(len(matrix)):
        for q in range(len(matrix[0])):
            if(matrix[j][q] != matrix[q][j]):
                flag1 = False
    for i in range (len(eigen_values)):
        if(eigen_values[i] < 0.0 or not(flag1)):
            flag = False

    return flag






def Cholesky_Decomposition(matrix,R,b,percision):
        print("\n Cholesky Method Steps : ")
        flag="correct"
        unique, infinite = Solution_test.IsUnique(matrix, b)
        if (not unique and not infinite):
            x = []
            flag = "System is INCONSISTENT"
            x.append(flag)
            return x,None
        if (infinite):
            x = []
            flag = "System has Infinite Solution"
            x.append(flag)
            return x,None

        t1=time.time()
        if(Hermitian_Positive_Semidefinite(matrix,R)):
    
            lower = [[0 for x in range(R)]for y in range(R)]

        # Decomposing a matrix
        # into Lower Triangular
            for i in range(R):
                for j in range(i + 1):
                    sum1 = 0

                    # summation for diagonals
                    #float('%.*g' % (precision, value) )
                    if (j == i):
                        for k in range(j):
                            sum1 += float('%.*g'%(percision,pow(lower[j][k], 2)))
                        lower[j][j] = float('%.*g'%(percision,math.sqrt(matrix[j][j] - sum1)))
                    else:

                        # Evaluating L(i, j)
                        # using L(j, j)
                        for k in range(j):
                            sum1 += float(('%.*g'%(percision,(lower[i][k] * lower[j][k]))))
                        if (lower[j][j] > 0):
                            lower[i][j] = float('%.*g'%(percision,(matrix[i][j] - sum1) /lower[j][j]))

            # Displaying Lower Triangular
            # and its Transpose
            # Ta2rib el matrix kolo
            #lower_rounded = np.round_(lower,decimals=percision)
            print("Lower Triangular\t\tTranspose")
            for i in range(R):

                # Lower Triangular
                transpose = [[0 for _ in range(R)]for _ in range(R)]
                for j in range(R):
                    print(lower[i][j], end="\t")
                print("", end="\t")

                # Transpose of
                # Lower Triangular
                for j in range(R):
                    transpose[i][j] = lower[j][i]
                    #transpose_rounded = np.round_(transpose,decimals=percision)
                    print(lower[j][i], end="\t")
                print("")
            t2 = time.time()
            x=solve(lower,transpose,b,percision)
            x.append(flag)
            return x,t2-t1
        else:
            print("Matrix is not a positive semidefinite matrix")
            flag="Matrix is not a positive semidefinite matrix"
            x=[]
            x.append(flag)
            return x,None


def solve(L,U,b,p):
    L = np.array(L)
    U = np.array(L.transpose())
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
    y[0] = float('%.*g'%(p,b[0]/L[0][0]))
    print("Y=",end="")
    print(y)
    print("\n")
    
    for i in range(1,len(y),1):
        sum0 = 0
        j=0
        for j in range(0,i,1):
            sum0 += float('%.*g'%(p,(L[i][j]*y[j])))
            y[i] = float('%.*g'%(p,((b[i] - sum0)/L[i][i])))
            #y = np.round_(y,decimals=p)
            
                
            print("Y= ",end="")
            print(y)

        
        
        
        
       #i = 1 - sum0 = l[1][0]*y[0]
        #y[i] = b[i] - sum / l[i][i]
        #i = 2 - sum0 = l[2][1]*y[1] + l[]
        

    print("Y= ",end="")    
    print(y)
    print("\n")
    x = [0 for _ in range(len(y))]
    
    x[len(y)-1] = float('%.*g'%(p,(y[len(y)-1]/U[len(y)-1][len(y)-1])))
    print("X= ",end="")
    print(x)
    for j in range(len(y)-2,-1,-1):
        sum1 = 0
        for q in range(j+1,len(y),1):
            sum1 += float('%.*g'%(p,(U[j][q]*x[q])))
            x[j] = float('%.*g'%(p,(y[j] - sum1)/U[j][j]))
            #x = np.round_(x,decimals=p)
            print("X= ",end="")
            print(x)
    print("X= ",end="")
    print(x)
    return x




    

