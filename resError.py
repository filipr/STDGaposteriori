from dolfin import *
import numpy as np 

from oneStep import oneStep

# compute the residual error 

# number of files 
N=0
error = []  

firstFile = '/home/filip/adgfem/RTN/export/nLevels.txt' 
with open(firstFile) as lines: 
      N, spaceLevels, timeLevels = [int(x) for x in next(lines).split()] # read first line 

print 'timeSteps = ' , N 
print 'spaceLevels = ', spaceLevels 
print 'timeLevels = ', timeLevels


# name of the files  
for i in range(N):
   if i<9: 
      fileName = '/home/filip/adgfem/RTN/export/rtn_00' + str(i+1)  + '.txt' 
   elif (i<99): 
      fileName = '/home/filip/adgfem/RTN/export/rtn_0' + str(i+1)  + '.txt' 
   else: 
      fileName = '/home/filip/adgfem/RTN/export/rtn_' + str(i+1)  + '.txt' 
   
   locErrorSquared = oneStep(fileName) 
   error.append( locErrorSquared ) 
   totError = np.sqrt( sum( error ) ) 
   print 'Local errors(',i+1,'):',  error
   print 'Total error:', totError 
