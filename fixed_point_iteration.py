import re
from re import sub
import time
from sympy import diff
import math
import matplotlib.pyplot as plt
import numpy as np

def substitute(equation,subs_value):
    equation = sub('[x]',str(subs_value), equation)
    return equation

def adjust_expression(exp):
    exp = exp.lower()    
    exp = exp.replace('^', '**')
    exp = exp.replace('log', 'math.log10')
    exp = exp.replace('sin', 'math.sin')
    exp = exp.replace('cos', 'math.cos')
    exp = exp.replace('tan', 'math.tan')
    exp = exp.replace('sinh', 'math.sinh')
    exp = exp.replace('cosh', 'math.cosh')
    exp = exp.replace('tanh', 'math.tanh')
    exp = exp.replace('pi', 'math.pi')
    exp = exp.replace('e', 'math.e')
    exp = exp.replace('sqrt', 'math.sqrt')
    
    return exp

def evaluate(expression,subs_value,SF):
    result=eval(substitute(adjust_expression(expression),subs_value))
    return float('%.*g' % (SF,result))

def abs_error(old,new,SF):
    error=abs(new-old) / abs(new)
    return(float('%.*g' % (SF,error)))

def fixedPointIteration(g,Xi, eps=0.00001, max_iter=50,SF=5,converge=True):

    t1=time.time()
    print("\033[1;32m"+"iter".ljust(10)+"Xi".ljust(20)+"Xi+1".ljust(20)+"E_apprx\x1B[0m")
    Xip=evaluate(g,Xi,SF)
    iters = 0
    error = abs_error(Xi,Xip,SF)
    div_count=0
    line=str(iters).ljust(10)+str(Xi).ljust(20)+str(Xip).ljust(20)+str(error)
    print(line)
    difference=abs(Xip-Xi)

    while error >= eps and max_iter > iters-1:
        Xi=Xip
        Xip=evaluate(g,Xi,SF)
        new_difference=abs(Xip-Xi)
        iters += 1
        prv_error=error
        error = abs_error(Xi,Xip,SF)
        line=str(iters).ljust(10)+str(Xi).ljust(20)+str(Xip).ljust(20)+str(error)
        print(line)
        
        if (new_difference>difference and not converge):
            div_count+=1
        else:
            div_count=0
        if div_count>2:
            t2=time.time()
            duration=t2-t1
            print("solution diverges")
            plot_function(g)
            return Xip,iters,error,"diverges",duration
        difference=new_difference
        

    t2=time.time()
    duration=t2-t1
    
    print("\033[1;94m"+"root",Xip)
    print("solution converges"+"\x1B[0m")
    plot_function(g,Xip-8,Xip+8)

    return Xip,iters,error,"converges",duration



def plot_function(g_x,lw=-10,up=10) :
    x = np.linspace(lw, up, 101)
    g_x=adjust_expression(g_x)
    y1=eval(g_x)
    y2 = x
    plt.title("fixed point iteration method")
    plt.axhline(y=0, color='k')
    plt.axvline(x=0, color='k')
    plt.plot(x,y1 ,'r-.')
    plt.plot(x, y2)
    plt.show()
   
def check_convergence_at_a_certain_root(equation,root):
    
    der=str(diff(adjust_expression(equation)))
    if(abs(eval(substitute(adjust_expression(der),root))) <1) :
         return True
    else :
         return False

if __name__ == '__main__':
    
    
    ##########testcases#############
    # expression="(2*x+3)^(1/2)"
    # expression="2*x*(1-x)"
    # initial_guess=0.2
    expression="-45/(x^3-18*x)"
    initial_guess=1.5
    # expression="6/(0.958*x^2-5.9*x+10.9)"
    # initial_guess=3.37
    # expression="(3.993*10^-4)/(x^2-0.165*x)"
    # expression="3/(2*x+3)"
    # expression="3/(x-2)"
    # expression="(x^2-3)/2"
    # initial_guess=4
    # initial_guess=0.05
    SFig=6
    eps=0.00001
    max_iter=50
    converge=check_convergence_at_a_certain_root(expression,initial_guess)

    if (converge):
        print("check:converges")
        
    else:
        print("check:diverges")
    root,iters,error,sol,t=fixedPointIteration(expression,initial_guess,eps,max_iter,SFig,converge)
   

    
    
    
    