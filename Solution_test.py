import numpy as np
def IsUnique(a,b):
   arank=np.linalg.matrix_rank(a)
   augmented=np.column_stack((a,b))
   augrank=np.linalg.matrix_rank(augmented)
   numberofVariables=len(a[0])
   unique=False
   infinite=False
   if(arank==augrank==numberofVariables):
      unique=True
   elif(arank==augrank and arank<numberofVariables):
      infinite=True
   return  unique,infinite