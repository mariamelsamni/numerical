import numpy as np
import time
import Solution_test
free_terms=['t','w','x','a','b','c','d','e','f','g','h','i','k','l','m','n','o','p','q','r','s','u','v','z','t1','w1','x1','a1','b1','c1','d1','e1','f1','g1','h1','i1','k1','l1','m1','n1','o1','p1','q1','r1','s1','u1','v1','z1']
free_terms_index=[]
Ans=[]
Ans.clear()
def decomposition(A,tolerance,n,b,error,scaling=False,precision=5):
    L=np.eye(n,dtype=float)
    # A=np.array(A)
    # A=np.round(A,precision)
    for i in range (n):
        for j in  range(n):
            A[i,j]=float('%.*g' % (precision, A[i,j]) )
    U=np.array(A,dtype=float)
    scale_factors=[1 for i in range(n)] 
    
    if (scaling):
        for i in range (n):
            scale_factors[i]= abs(A[i,0])                 
            for j in range (1,n):
                if ( abs(A[i,j]) > scale_factors[i] ):    
                    scale_factors[i]= abs(A[i,j])
    print("scale factors",scale_factors)
    
    print("\n\033[1;32m=> forward elimination to decompose coefficient matrix into upper and lower "+",\x1B[0m")
    for c_r in range (n):
        A,b=pivoting(A=A,c_r=c_r,tolerance=tolerance,scale_factors=scale_factors,b=b,n=n,L=L,U=U)
        if (A[c_r,c_r]==0):#if the pivot after applying partial pivoting is zero  
            print("scale factor ",scale_factors[c_r])          
            free_terms_index.append(c_r)
            error=2
            print("infinite",c_r)
            continue
        
        if ((abs(A[c_r,c_r])/scale_factors[c_r])<tolerance):
            print("error near singular")
            return A,b,1

        for i in range (c_r+1,n):
            error=1
            factor=float('%.*g' % (precision, A[i,c_r]/A[c_r,c_r]) )
            A[i,c_r]=factor
                        
            print("\x1B[4m"+"row"+str(i)+" = row",str(i)+" - ",str(factor)+ " * row", str(c_r)+",\x1B[0m" )
            for j in range(c_r+1,n):
                A[i,j]-=float('%.*g' % (precision, factor*(A[c_r,j]) ))
                A[i,j]=float('%.*g' % (precision, (A[i,j]) ))
            L,U=get_L(A,i,L,U)      
            print("\nU \n",U,"\nL \n",L)
    if (error!=2):
        if ((abs(A[n-1,n-1])/scale_factors[n-1])<tolerance):
                print("error near singular")
                return A,b,1

    L,U=get_L(A,i,L,U)
    print("\nfinal U \n",U,"\nfinal L \n",L)
    print("A\n",A)
    return A,b,0


def substitution(A,b,tolerance,n,precision=5):
    A=np.array(A,dtype=float)
    print("\nb\n",b)
    print("\n\033[1;32m=> forward substitution to solve Lz=b "+",\x1B[0m")
    z=np.array([0 for i in range (n)],dtype=float)
    x=([0 for i in range (n)])    
    for i in range (n):
        b[i]=float('%.*g' % (precision, (b[i]) ))
    print ("z[0] = "+ str(b[0]))
    z[0]=b[0]
    for row in range (1,n):
        sum=b[row]
        print_='z['+str(row)+"]="+str(sum)
        for col in range(row):
            sum-=float('%.*g' % (precision, (A[row,col]*z[col]) ))
            sum=float('%.*g' % (precision, (sum) ))
            print_=print_+"-"+str(A[row,col])+"*"+str(z[col])
        print(print_)
        z[row]=sum
    print(">>>z",z)
    print("\n\033[1;32m=> backward substitution to solve Ux=z "+",\x1B[0m")
    print ("x["+str(n-1)+"] = "+ str(z[n-1])+"/"+str(A[n-1,n-1]))
    print("\nz[n-1]\n",z[n-1])
    if (A[n-1][n-1]==0  ):
        free_term_not_yet=True
        x[n-1]=free_terms[n-1]
        
    else:
        x[n-1]=float('%.*g' % (precision, (z[n-1]/A[n-1,n-1]) ))
        free_term_not_yet=False
    for row in range(n-2,-1,-1):
        
        if (row in free_terms_index):
            free_term_not_yet=True
            x[row]=free_terms[row]
            continue
        if (not free_term_not_yet):
            sum=z[row]
            print_='x['+str(row)+"]=("+str(sum)
            for col in range (row+1,n):
                sum-=float('%.*g' % (precision, (A[row,col]*x[col]) ))
                sum=float('%.*g' % (precision, (sum) ))
                print_=print_+"-"+str(A[row,col])+"*"+str(x[col])
            x[row]=float('%.*g' % (precision, sum/A[row,row] ))
            print_=print_+")/"+str(A[row,row])
            print(print_)
        else:
            sum=float('%.*g' % (precision, z[row]/A[row,row] ))
            print_='x['+str(row)+"]=("+str(sum)
            attached_str=""
            for col in range (row+1,n):
                if (isinstance(x[col],str)): 
                # if (x[col]in free_terms):                    
                    attached_str+="-"+str(float('%.*g' % (precision,((A[row,col])/A[row,row]))))+x[col]
                    continue
                print("elmoshekla\n",A[row,col])
                print("elmoshekla\n",x[col])
                sum-=float('%.*g' % (precision,float('%.*g' % (precision, (A[row,col])*x[col]) ))/A[row,row])
                print_=print_+"-"+str(A[row,col])+"*"+str(x[col])
            x[row]=str(sum)+attached_str
            print_=print_+attached_str+")"
            print(print_)
    for i in range(len(x)):
        Ans.append(x[i])
    print("x", x)
    print("\033[1;94m"+">>>x "+str(x)+"\x1B[0m")
        
def get_L(A,c_r,L,U):
    
    U=np.array(A)
    for i in range (c_r+1):
        for j in range (i):
                U[i,j]=0
                L[i,j]=A[i,j]
    
    return L,U


def pivoting(A,c_r,tolerance,scale_factors,b,n,L,U):

    if (scale_factors[c_r]==0):
        pivot=0
    else:
        pivot=abs(A[c_r,c_r])/scale_factors[c_r]
    row_of_pivot=0                  
    for i in range (c_r,n):
        if (scale_factors[i]==0):
            continue
        if ( abs(A[i,c_r])/scale_factors[i] > pivot):
            pivot= abs(A[i,c_r])/scale_factors[i]
            row_of_pivot=i
    if (row_of_pivot):
        print("swap row ", c_r," with row  ", row_of_pivot,"in both matrices and the result (b) matrix")
        A[[c_r,row_of_pivot]]=A[[row_of_pivot,c_r]]
        b[[c_r,row_of_pivot]]=b[[row_of_pivot,c_r]]
        L,U=get_L(A,row_of_pivot,L,U)
        print("\nU \n",U,"\nL \n",L)
        print("\nb\n",b)
        print("-------------------------------------------------")

    return A,b
            
def LU_doolittle(A,b,tolerance,n,precision=3,scale=False):
    print("\n Doolittle Method Steps : ")
    Ans.clear()
    unique, infinite = Solution_test.IsUnique(A, b)
    if (not unique and not infinite):
        x = []
        flag = "System is INCONSISTENT"

        return x,flag,None
    t1=time.time()
    x=np.array([0 for i in range (n)])
    error="correct"
    A=np.array(A,dtype=float)
    b=np.array(b,dtype=float)
    # isZero = np.all(A == 0, axis=1)
    # count=0
    # for i in range (n):
    #     if (isZero[i] ):
    #         count+=1
    # if (count==1):
    #     print("infinite")
    #     t2=time.time()
    #     t=t2-t1
    #     flag="infinite"
    #     return Ans,flag,t
    A,b,error=decomposition(A=A,tolerance=tolerance,n=n,b=b,error=error,scaling=scale,precision=precision)
    
    if error==1:
        print("system can't be solved using LU doolittle method (near singular)")
        error="system can't be solved using LU doolittle method (near singular)"
    else:
        x=substitution(A=A,b=b,tolerance=tolerance,n=n,precision=precision)
        error="correct"
    t2=time.time()
    t=t2-t1
    print("time: ",t)
    return(Ans,error,t)
    

if __name__ == '__main__':

    n=5  #the number of equations
    # coefficients=[[2,1,5],[4,4,-4],[0,0,0]]
    # equationsResults=np.array([5,0,0],dtype=float)
   
    # coefficients=[[0,0],[0,0]]
    # equationsResults=np.array([0,0],dtype=float)
    coefficients=[[1,3,2,4,-3],[2,6,0,-1,-2],[0,0,6,2,-1],[1,3,-1,4,2],[0,0,0,0,0]]
    equationsResults=np.array([-7,0,12,-6,0],dtype=float)
    #>>>>>>>>>>>>>>>>>>>>TEST CASES<<<<<<<<<<<<<<<<<<<<<<<<<    
    # coefficients=[[2,3,5,8],[2,3,6,9],[2,3,7,10],[2,3,7,11]]
    # equationsResults=np.array([18,20,22,23],dtype=float)
    # coefficients=[[2,3,5,8],[2,3,5,8],[2,3,5,8],[2,3,5,8]]
    # equationsResults=np.array([18,18,18,18],dtype=float)
    # coefficients=[[0,2],[3,0]]
    # equationsResults=np.array([2,3],dtype=float)
    # coefficients=[[1,-1],[5,-5]]
    # equationsResults=np.array([3,15],dtype=float)
    # coefficients=[[0,0],[0,0]]
    # equationsResults=np.array([3,15],dtype=float)
    # coefficients=[[0,1],[0,2]]
    # equationsResults=np.array([1,2],dtype=float)
    # coefficients=[[1,2],[0,0]]
    # equationsResults=np.array([0,0],dtype=float)
    # coefficients=[[6,15,55],[15,55,225],[55,225,979]]
    # equationsResults=np.array([5,0,6],dtype=float)
    # coefficients=[[25,5,1],[64,8,1],[144,12,1]]
    # equationsResults=np.array([1,-2.56,3.2],dtype=float)
    # coefficients=[[12,3,-5],[1,5,3],[3,7,13]]
    # equationsResults=np.array([1,28,76],dtype=float)
    # coefficients=[[4,2,1],[-1,2,0],[2,1,4]]
    # equationsResults=np.array([11,3,16],dtype=float)
    # coefficients=[[1,2],[2,1]]
    # equationsResults=np.array([3,3],dtype=float)
    # coefficients=[[1.2,1,1,1,1],[1,2,1,1,1],[1,1,2,1,1],[1,1,1,2,1],[1,1,1,1,2]]
    # equationsResults=np.array([4,5,6,7,8],dtype=float)
    #>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    tolerance=0.0000005
    precision=5

    error,x,t=LU_doolittle(A=coefficients,b=equationsResults,tolerance=tolerance,n=n,precision=precision,scale=True)
