import time
import functools
import numpy as np
import collections
def isZeroOnDiagonal(a):#check if the matrix has zero on diagonal
    zeroOnDiagonal=False
    for i in range(len(a)):
        if a[i][i]==0 :
           zeroOnDiagonal=True
           break
    return zeroOnDiagonal
def AvoidZeroOnDiagonalIfPossible(a,b): #Re-arrange the matrix to avoid zeros on the diagonal if possible
    if (not isZeroOnDiagonal(a)): return
    for i in range(len(a)):
        if a[i][i]!=0:continue
        for j in range(len(a)):
            if a[j][i]!=0: # if the element in row j (in the same column of the diagonal element[i][i] ) doesn't equal zero then switch rows i and j
                t=a[i]
                a[i]=a[j]
                a[j]=t
                t = b[i]  #switch the equations results matrix rows  as well
                b[i] = b[j]
                b[j] = t

        if(not isZeroOnDiagonal(a)): return
def isDiagonalElementMaxInRow(a,i): #check if the absolute value of the diagonal element in a specific row is larger than the sum of the absolute value of the others elements of the row
    greater=False
    greater_equal=False
    absoluteSum=0
    for j in range(len(a[i])):
        absoluteSum=absoluteSum+abs(a[i][j])
    if abs(a[i][i])>=absoluteSum-abs(a[i][i]):
        greater_equal=True
        if abs(a[i][i])>absoluteSum-abs(a[i][i]):
            greater=True

    return greater_equal,greater
def isDiagonalDominant(a): #check if the matrix is diagonally dominant
    DD=True
    G=False
    for i in range(len(a)):
        greater_equal,greater=isDiagonalElementMaxInRow(a,i)
        if(not greater_equal):
            DD=False
        elif(greater):
            G=True

    return (DD and G)
def isMaxElementInRowLargerThanSumOfTheRest(a,i): #check if there is an element in the row whose absolute value is the larger than the sum of the absolute value of the others elements in the row
    absoluteSum = 0
    for j in range(len(a[i])):
        absoluteSum = absoluteSum + abs(a[i][j])
    if (max(a[i])>=abs(min(a[i]))):# handling the absolute value case
        if max((a[i]))>=absoluteSum-max((a[i])):
            return True
        else:return False
    if(abs(min(a[i]))>=max(a[i])):
        if abs(min(a[i]))>=absoluteSum-abs(min(a[i])):
            return True
        else:return False
def indexOfMax(a,i): #return index of the largest element in row whose absolute value is the larger than the sum of the absolute value of the others elements in the row
    if (max(a[i]) >= abs(min(a[i]))):
        return a[i].index(max(a[i]))
    if (abs(min(a[i])) >= max(a[i])):
        return a[i].index(min(a[i]))
def rearrangeMatrix(a,b):#re-arrange the matrix to transform it into diagonally dominant if possible
    counter=0
    while counter<len(a):
       if (isDiagonalDominant(a)): return
       for i in range(len(a)):
        greater_equal,greater=isDiagonalElementMaxInRow(a,i)
        if(greater_equal):continue
        elif(isMaxElementInRowLargerThanSumOfTheRest(a,i)):
                j=indexOfMax(a,i)
           #if the row of max element already doesn't have  a diagonally dominant element then switch rows
                t=a[i]
                a[i]=a[j]
                a[j]=t
                t=b[i]
                b[i]=b[j]
                b[j]=t
       counter=counter+1

    AvoidZeroOnDiagonalIfPossible(a,b)

def Calculate_X_Seidel(x,a,b,i,precision):
    x[i]=float('%.*g' % (precision,b[i]))
    for j in range(len(a)):
        if j!=i:
            # float('%.*g' % (precision,))
            # x[i]=np.round(x[i],precision)-np.round(a[i][j],precision)*np.round(x[j],precision)
            x[i] =float('%.*g' % (precision,float('%.*g' % (precision,x[i])) - float('%.*g' % (precision,a[i][j])) * float('%.*g' % (precision,x[j]))))
    # x[i]=np.round((x[i]/a[i][i]),precision)
    x[i] = float('%.*g' % (precision,(x[i]/float('%.*g' % (precision,a[i][i])))))
    if (x[i] == -0):
        x[i] = abs(x[i])
def Calculate_X_Jacobi(x_new,x_old,a,b,i,precision):
    x_new[i]=b[i]
    for j in range(len(a)):
        if j!=i:
            # x_new[i]=np.round(x_new[i],precision)-np.round(a[i][j],precision)*x_old[j]
            x_new[i] = float('%.*g' % (precision,float('%.*g' % (precision, x_new[i])) - float('%.*g' % (precision, a[i][j])) * float('%.*g' % (precision, x_old[j]))))
    # x_new[i]=np.round(x_new[i]/a[i][i],precision)
    x_new[i] = float('%.*g' % (precision, (x_new[i] / float('%.*g' % (precision, a[i][i])))))
    if(x_new[i]==-0):
        x_new[i]=abs(x_new[i])
def Seidel(a,b,guess,precision,rearrange_condition,tolerance=-1,NumberofItirations=100):
    ErrorsList = []
    ErrorsList.clear()
    print("\n Seidel Method Steps : ")
    x = guess.copy()
    x_old=guess.copy()
    flag = "correct"
    if (rearrange_condition == "Re-arrange Matrix"):
        rearrangeMatrix(a, b)
        print("Matrix after re-arrangement :")
        for i in range(len(a)):
         print(a[i])
    if isZeroOnDiagonal(a):
        print("Cannot Solve : Diagonal Contains ZERO")
        flag="Cannot Solve : Diagonal Contains ZERO"
        x.append(flag)
        return x,None
    tolerance_meeted_counter=0
    t1 = time.time()
    Itirations=0
    while Itirations<100:
          if Itirations==NumberofItirations:
              break
          print("itiration ",Itirations+1," :",)
          for i in range(len(x)):
              x_prev=x[i]
              Calculate_X_Seidel(x,a,b,i,precision)
              if(x[i]!=0):
               error=float('%.*g' % (precision,abs(float('%.*g' % (precision,(x[i]-x_prev)))/x[i])))
               if error<=tolerance:
                  tolerance_meeted_counter=tolerance_meeted_counter+1
                  ErrorsList.append(error)
          if x_old == x:
              print("x=", x)
              if (tolerance == -1): tolerance = 0
              print("tolerance", tolerance, "meeted")
              print("Absolute Relative Error = 0")
              break
          if tolerance_meeted_counter==len(x):
              print("x=",x)
              print("tolerance", tolerance, "meeted")
              print("Absolute Relative Error = ", max(ErrorsList))
              break
          else:
              tolerance_meeted_counter=0
              ErrorsList.clear()
          Itirations=Itirations+1
          print("x= ",x)
          x_old=x.copy()
    t2 = time.time()
    if Itirations==100:
        print("System didn't converge")
        flag="System didn't converge"
        x.append(flag)
        return x,None
    else:
     x.append(flag)
     return x,t2-t1
def Jacobi(a,b,guess,precision,rearrange_condition,tolerance=-1,NumberofItirations=100):
    ErrorsList=[]
    ErrorsList.clear()
    print("\n Jacobi Method Steps : ")
    x_new = [0 for i in range(len(a[0]))]
    x_old=guess.copy()
    flag="correct"
    if(rearrange_condition=="Re-arrange Matrix"):
     rearrangeMatrix(a,b)
     print("Matrix after re-arrangement :")
     for i in range(len(a)):
        print(a[i])
    if isZeroOnDiagonal(a):
        print("Cannot Solve : Diagonal Contains ZERO")
        flag="Cannot Solve : Diagonal Contains ZERO"
        x_new.append(flag)
        return x_new,None
    tolerance_meeted_counter=0
    t1 = time.time()
    Itirations=0
    while Itirations<100:
          if Itirations == NumberofItirations:
            break
          print("itiration ",Itirations+1," :",)
          for i in range(len(x_new)):

            Calculate_X_Jacobi(x_new,x_old,a,b,i,precision)
            if(x_new[i]!=0):
              # error=abs((x_new[i]-x_old[i])/x_new[i])
              error = float('%.*g' % (precision, abs(float('%.*g' % (precision, (x_new[i] - x_old[i]))) / x_new[i])))
              if error<=tolerance:
                  tolerance_meeted_counter=tolerance_meeted_counter+1
                  ErrorsList.append(error)
          if x_old == x_new:
              print("x_old = ", x_old)
              print("x_new = ", x_new)
              if(tolerance==-1):tolerance=0
              print("tolerance", tolerance, "meeted")
              print("Absolute Relative Error = 0")
              break
          if tolerance_meeted_counter==len(x_new):
              print("x_old = ", x_old)
              print("x_new = ", x_new)
              print("tolerance",tolerance ,"meeted")
              print("Absolute Relative Error = ",max(ErrorsList))
              break
          else:
              tolerance_meeted_counter=0
              ErrorsList.clear()
          Itirations=Itirations+1
          print("x_old = ",x_old)
          print("x_new = ",x_new)
          x_old=x_new.copy()
    t2 = time.time()

    if Itirations==100:
        print("System didn't converge")
        flag="System didn't converge"
        x_new.append(flag)
        return x_new,None
    else:
         x_new.append(flag)
         return x_new,t2-t1



