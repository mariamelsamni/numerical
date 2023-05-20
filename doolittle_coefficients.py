import numpy as np
free_terms=['t','w','x','a','b','c','d','e','f','g','h','i','k','l','m','n','o','p','q','r','s','u','v','z']
free_terms_index=[]
Ans=[]
Ans.clear()
import time
def clone(A,n):
    U=[[ 0 for i in range (n)] for i in range (n)]
    for i in range(n):
        for j in range(n):
            U[i][j]=A[i][j]
    return U


def decomposition(A,n,b,error):
    L=np.eye(n)
    L=L.tolist()
    U=clone(A,n)
    
    
    
    print("\033[1;32m=> forward elimination to decompose coefficient matrix into upper and lower "+",\x1B[0m")
    for c_r in range (n):
        
        for i in range (c_r+1,n):
           
            if (A[i][c_r]==A[c_r][c_r]):
                print ("infinite or no solution")
                return A,b,1
            factor=str(A[i][c_r])+"/"+str(A[c_r][c_r])
            A[i][c_r]=factor
            print("\x1B[4m"+"row"+str(i)+" = row",str(i)+" - ",str(factor)+ " * row", str(c_r)+",\x1B[0m" )
            for j in range(c_r+1,n):
                A[i][j]+="-"+factor+ "*" +str(A[c_r][j])
                U[i][j]+="-"+factor+ "*" +str(A[c_r][j])
            L,U=get_L(A,i,L,U,n)
            print("\nU \n",U,"\nL \n",L)
          
    return A,b,error


def substitution(A,b,n):
    print("\n\033[1;32m=> forward substitution to solve Lz=b "+",\x1B[0m")
    
    z=[0 for i in range (n)]
    x=([0 for i in range (n)])
    print ("z[0] = "+ str(b[0]))
    z[0]=b[0]
    
    for row in range (1,n):
        sum=b[row]
        print_1='z['+str(row)+"]="
        print_2=str(sum)
        for col in range(row):
            print_2=print_2+"-"+str(A[row][col])+"*"+str(z[col])
        print(print_1+print_2)
        z[row]=print_2

    print("z",z)
    print("\n\033[1;32m=> backward substitution to solve Ux=z "+",\x1B[0m")
    print ("x["+str(n-1)+"] = "+ str(z[n-1])+"/"+str(A[n-1][n-1]))
   
    
    x[n-1]=str(z[n-1])+"/"+str(A[n-1][n-1])
    
    
    for row in range(n-2,-1,-1):
        
        sum=z[row]
        print_2=z[row]
        print_='x['+str(row)+"]="
        for col in range (row+1,n):
            print_2=print_2+"-"+str(A[row][col])+"*"+str(x[col])
        x[row]=print_2
        print_=print_+")/"+str(sum)+"/"+str(A[row][row])
        print(print_+print_2)
    for i in range(len(x)):
        Ans.append(x[i])
    print("x",x)
        
def get_L(A,c_r,L,U,n):

    U=clone(A,n)
    for i in range (c_r+1):
        for j in range (i):
                L[i][j]=A[i][j]
                U[i][j]=0

    return L,U

            
def LU_doolittle(A,b,n):        
    Ans.clear()
    x=[0 for i in range (n)]
    error=0
    t1=time.time()
    A,b,error=decomposition(A=A,n=n,b=b,error=error)
    
    if error==1:
        print("system can't be solved using LU doolittle method ")
        error="system can't be solved using LU doolittle method "
    else:
        x=substitution(A=A,b=b,n=n)
        error="correct"
    t2=time.time()
    return(Ans,error,t2-t1)
    

