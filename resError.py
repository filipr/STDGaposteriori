from dolfin import *
import numpy as np 

from oneStep import readRHS, dualError
from tex import makeTex


directoriesCompute = [ \
               'case02_p02_q02_nelem000512_steps00020' \
                ]
                
directoriesTex = ['case02_p01_q01_nelem000128_steps00005', \
               'case02_p01_q01_nelem000512_steps00010', \
               'case02_p01_q01_nelem002048_steps00020', \
               'case02_p02_q01_nelem000032_steps00002', \
               'case02_p02_q01_nelem000128_steps00010', \
               'case02_p02_q01_nelem000512_steps00020', \
               'case02_p02_q02_nelem000512_steps00020' \
                ]


# compute the dual error of the residual for all directories
# then save the results to files with the same name as the directory 
def compute(directories):
   for i in range(len(directories)): 
      line, lineEtas, case, eps, ipg = dualError( directories[i] ) 
      
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
compute(directoriesCompute) 

# Make tex table from the given files !!!
makeTex(directoriesTex)
   
   







