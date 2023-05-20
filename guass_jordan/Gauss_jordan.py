
import numpy as np
import sys
from copy import deepcopy
import time
import Solution_test




# guass jordan for numbers
def guass_jordan_numbers(a, b, scale, precision):
    print("\n Guass Jordan Method Steps : ")
    unique, infinite = Solution_test.IsUnique(a, b)
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
    start = time.time()
    n = len(b)
    # xx = len(a)
    # yy = len(a[0])
    iterations = 0
    # if xx != yy:
    #     sys.exit("Sorry cant proceed,Number of equations should be equal number of variables")
    for i in range(n):
        if (scale):
            # scaling

            sama = [[0] * n] * n
            sama = deepcopy(a)
            for st in range(i+1,n):
                it = a[st][0]
                for sz in range(n):
                    if it < a[st][sz]:
                        it = a[st][sz]
                if it < b[st]:
                    it = b[st]
                for sa in range(n):
                  if(it!=0):
                    sama[st][sa] = sama[st][sa] / it

            # partial pivoting each loop of rows
            flag = i
            value = np.fabs(sama[i][i])
            if (i != (n - 1)):
                for j in range(n):
                    if value < np.fabs(sama[j][i]):
                        flag = j
                        value = np.fabs(sama[j][i])
            print("R", i, "↔", "R", flag)
            for k in range(n):
                temp = a[i][k]
                a[i][k] = a[flag][k]
                a[flag][k] = temp
            temp2 = b[flag]
            b[flag] = b[i]
            b[i] = temp2
        else:
            flag = i
            value = np.fabs(a[i][i])
            if (i != (n - 1)):
                for j in range(n):
                    if value < np.fabs(a[j][i]):
                        flag = j
                        value = np.fabs(a[j][i])
            print("R", i, "↔", "R", flag)
            for k in range(n):
                temp = a[i][k]
                a[i][k] = a[flag][k]
                a[flag][k] = temp
            temp2 = b[flag]
            b[flag] = b[i]
            b[i] = temp2

        # division of pivot(pivot=1)
        pivot = a[i][i]
        print("R", i, "↔", "R", i, "/", pivot)
        b[i] = float('%.*g' % (precision, float('%.*g' % (precision, b[i]) )/float('%.*g' % (precision,pivot) )) )
        for t in range(n):
            a[i][t] = float('%.*g' % (precision, float('%.*g' % (precision, a[i][t]) )/float('%.*g' % (precision,pivot) )) )
        # elimination to echelon form
        for s in range(n):
            if s == i or a[s][i] == 0:
                continue
            ratio = float('%.*g' % (precision,a[s][i]) )
            b[s] = float('%.*g' % (precision,float('%.*g' % (precision,b[s]))-float('%.*g' % (precision,ratio) )*float('%.*g' % (precision,b[i]) )) )
            for q in range(i, n):
                a[s][q] = float('%.*g' % (precision,float('%.*g' % (precision,a[s][q]))-float('%.*g' % (precision,ratio) )*float('%.*g' % (precision,a[i][q]) )) )
            print("R", s, "↔", "R", s, "-", ratio, "*", "R", i)

        # printing matrix after each step
        print("Matrix after", iterations + 1, "iterations:")
        iterations += 1
        for r in range(n):
            for y in range(n):
                print(a[r][y], end="      ")
            print()
        print()
    end = time.time()
    b.append("correct")
    return b, end - start


##########################################################################################################


# guass jordan for characters
def guass_jordan_characters(a, b):
    start2 = time.time()
    n = len(b)
    iterations = 0
    # xx = len(a)
    # yy = len(a[0])
    # if xx != yy:
    #     sys.exit("Sorry cant proceed,Number of equations should be equal number of variables")
    for i in range(n):
        # division of pivot(pivot=1)
        pivot = a[i][i]
        print("R", i, "↔", "R", i, "/", pivot)
        b[i] = '(' + b[i] + '/' + pivot + ')'
        for t in range(n):
            a[i][t] = '(' + a[i][t] + '/' + pivot + ')'
        # elimination to echelon form
        for s in range(n):
            if s == i or a[s][i] == 0:
                continue
            ratio = a[s][i]
            b[s] = '(' + b[s] + '-' + ratio + '*' + b[i] + ')'
            for q in range(i, n):
                a[s][q] = '(' + a[s][q] + '-' + ratio + '*' + a[i][q] + ')'
            print("R", s, "↔", "R", s, "-", ratio, "*", "R", i)

        # printing matrix after each step
        print("Matrix after", iterations + 1, "iterations:")
        iterations += 1
        for r in range(n):
            for y in range(n):
                print(a[r][y], end="      ")
            print()
        print()
    end2 = time.time()
    b.append("correct")
    return b, end2 - start2


#driver code
if __name__ == '__main__':
    a=[['a','b','c'],['d','e','f'],['g','h','i']]
    b=['q','r','s']
    # c=[[1,1,2],[3,4,7],[2,1,4]]
    c=[[1,1,-3,1],[-5,3,-4,1],[1,0,2,-1],[1,2,0,0]];
    # d=[8,10,6]
    d=[2,0,1,12]
    Y=guass_jordan_numbers(c,d,1,5)
    print("The answer is:")
    print(Y)
    print(" ")
    X=guass_jordan_characters(a,b)
    print("The answer is:")
    print(X)