import numpy as np
import time
Ans=[]
Ans.clear()
def pivoting(arr1,k,result,precision):
    #assume the first element is the max
    max = [k,float('%.*g' % (precision,abs(arr1[k][k])))]
    # compare the pivot with all the numbers under it
    for i in range(k+1,len(arr1)):
        if(float('%.*g' % (precision,abs(arr1[i][k])))>float('%.*g' % (precision,max[1]))):
            max[0] = i
            max[1] = float('%.*g' % (precision,abs(arr1[i][k])))
    # if there is a number greater than the pivot ... then change the two rows
    if(max[0]>k):
        print("pivoting by swap R"+str(k+1)+",","R"+str(max[0]+1))
        temp = arr1[k]
        arr1[k] = arr1[max[0]]
        arr1[max[0]] = temp
        temp = result[k]
        result[k] = result[max[0]]
        result[max[0]] =  temp
        print(arr1)
        print(result)
        print("_______________________________________________________________________________________________________________________")
    return
        
        
def backSubstitution(arr1,result,precision):
    flag="correct"
    print("Back Substitution")
    l = len(arr1)-1
    x = np.ones(len(result))
    # check that the diagonal elements are not zeros ... if one of them is zero ... then there are infinite number of solutions
    for i in reversed(range(len(arr1))):
        check = True
        if(arr1[i][i]==0):
            for j in range(i+1,len(arr1)):
                if(arr1[i][j]==0):
                    continue
                else:
                    check = False
                    x[i] = 0 #new
            if(check==True):
                if(result[i]==0):
                    x[i] = 0 #new
                else:
                    print("The System has no solution")
                    y=[]
                    flag="The System has no solution"
                    return y,flag


    if(x[l]!=0):
        # calculate the last x 
        x[l] = float('%.*g' % (precision,result[l] / arr1[l][l]))
    else:
        flag="infinite"
        print("Infinte number of solutions")
    print("X"+str(l+1),"=",x[l])
    #loop to calculate Xi where i = [n-1,1]
    for i in reversed(range(l)):
        sum = 0
        for j in range(i+1,len(arr1)):
            sum = float('%.*g' % (precision, sum +  float('%.*g' % (precision,arr1[i][j]*x[j]))))
        if(x[i]!=0):
            x[i] = float('%.*g' % (precision,(result[i]-sum)/arr1[i][i]))
        print("X"+str(i+1),"=",x[i])

    return x,flag

#function to get the maximum absolute number in the row
def maximum(arr):
    max = abs(arr[0])
    for i in range(len(arr)):
        if(abs(arr[i])>max):
            max = abs(arr[i])
    if(max==0):
        return 1
    else:
        return max


def pivotingWithScalling(arr1,k,result,precision):
    max = [k,float('%.*g' % (precision,abs(arr1[k][k]/maximum(arr1[k]))))]
    # compare the pivot with all the numbers under it
    for i in range(k+1,len(arr1)):
        if(float('%.*g' % (precision,abs(arr1[i][k]/maximum(arr1[i]))))>max[1]):
            max[0] = i
            max[1] = float('%.*g' % (precision,abs(arr1[i][k]/maximum(arr1[i]))))
    # if there is a number greater than the pivot ... then change the two rows
    if(max[0]>k):
        print("pivoting by swap R"+str(k+1)+",","R"+str(max[0]+1))
        temp = arr1[k]
        arr1[k] = arr1[max[0]]
        arr1[max[0]] = temp
        temp = result[k]
        result[k] = result[max[0]]
        result[max[0]] =  temp
        print(arr1)
        print(result)
        print("_______________________________________________________________________________________________________________________")
    return






def GuassElimination(arr1,result,scalling,precision):
    for i in range(len(arr1)):
        for j in range(len(arr1)):
            arr1[i][j] = float('%.*g' % (precision,arr1[i][j]))
    for i in range(len(result)):
        result[i] = float('%.*g' % (precision,result[i]))
    Ans.clear()
    print("\n Guass Elimination Method Steps : ")
    start = time.time()
    print(arr1)
    print(result)
    print("_______________________________________________________________________________________________________________________")
    for k in range (len(arr1)-1):
        # for each pivot make the pivoting function ... you can make pivoting with scalling and without
        if(scalling==True):
            pivotingWithScalling(arr1,k,result,precision)
        else:
            pivoting(arr1,k,result,precision)
        for i in range(k+1,len(arr1)):
            #check if the elemnt is not zero ... then make it zero
            if(arr1[i][k]==0): 
                continue
            # calculate the factor for the row
            factor = arr1[i][k]/arr1[k][k]
            # print the operation 
            print(-factor,"*","R"+str(k+1),"+ R"+str(i+1))
            #loop to make the operation
            for j in range(len(arr1)):
                arr1[i][j] = float('%.*g' % (precision,float('%.*g' % (precision, arr1[k][j]*(-factor))) + arr1[i][j] ))    
            result[i] = float('%.*g' % (precision,result[i] - factor*result[k]))
            print(arr1)
            print(result)
            print("_______________________________________________________________________________________________________________________")
    # making the substitution function after getting the pivots
    x,flag= backSubstitution(arr1,result,precision)
    end = time.time()
    print("Time to calculate =",end-start,"sec")
    print("_______________________________________________________________________________________________________________________")
    for i in range(len(x)):
        Ans.append(x[i])
    Ans.append(flag)
    return Ans,end-start




        








        










        




