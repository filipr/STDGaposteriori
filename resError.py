# from dolfin import *
# import numpy as np 

from oneStep import dualError
# from tex import makeTex

# compute the dual error of the residual for all directories
# then save the results to files with the same name as the directory 
def compute(directories, computeDual):
   for i in range(len(directories)): 
      line, lineEtas, case, eps, ipg = dualError( directories[i] , computeDual) 
      
      # write to file
      fileName = directories[i] 
      with open(fileName, 'w') as f:
            f.write( '{:5}{:2}{:5}{:2}{:5}{:2}'.format(case,'  ', eps,'  ', ipg, ' \n') )
            for x in line: 
               f.write( '{:}{:2}'.format(x, '  ') )
            f.write('\n')
            for x in lineEtas: 
               f.write( '{:}{:2}'.format(x, '  ') )
            f.close() 
   return  
   
   
#Compute the errors and save them to files 
# compute(directoriesCompute) 

# Make tex table from the given files !!!
#makeTex(directoriesTex)
   
   







