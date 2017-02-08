from dolfin import *
import numpy as np 

from oneStep import oneStep 
from tex import makeTex

def dualError(directory): 
   '''
   N =  number of files, i.e. number of time steps
   spaceLevels = number of little elements per coarse triangle 
   spaceLevels = number of little time steps pre one coarse time step 
   nelem = #elements in the coarse mesh 
   h = coarse H
   tau = coarse time step length 
   p = space pol. degree
   q = time pol. degree
   case = isca case in Adgfem 
   eps = epsilon in front of the diffusion member 
   ipg = SIPG,... 
   etaR, etaF, etaT, etaNC = error estimators (Rezidual, Flux, timeRadau, NonConformity)
   eta = discretizazion estimator (etaR + etaF + etaT) 
   etaApprox = approximate of the error upper bound? 
   '''

   folder = '/home/filip/adgfem/RTN/export/'+ directory + '/'
   firstFile = folder + 'initDualProblem.txt' 
   
   with open(firstFile) as lines: 
         # N 
         N, spaceLevels, timeLevels = [int(x) for x in next(lines).split()] # read first line 
         nelem, h , tau, p, q    = [float(x) for x in next(lines).split()]
         case, eps, ipg          = [x for x in next(lines).split()] 
         etaR, etaF, etaT, etaNC = [float(x) for x in next(lines).split()]
         eta , errApprox         = [float(x) for x in next(lines).split()]
         lines.close() 
   
   print 'timeSteps = ' , N 
   print 'spaceLevels = ', spaceLevels 
   print 'timeLevels = ', timeLevels 
   error = []
   # name of the files  
   for i in range(N): 
      fileName = folder + 'rtn_' + str(i+1).zfill(4) + '.txt'
      error.append( oneStep(fileName) )  
      print 'Total error (step', i,'): ', np.sqrt( sum( error ) )
   totError = np.sqrt( sum( error ) ) 
   print 'Total error:', totError 
   
   iEff = eta / totError 
   iEffTot = ( eta + etaNC) / ( totError + etaNC )   
   iEffPseudo = eta / errApprox 
   
   line = [nelem, h, tau, p, q, totError, eta, errApprox, etaNC, iEff, iEffTot, iEffPseudo ]
   lineEtas = [etaR, etaF, etaT, etaNC] 
   
   return line, lineEtas, case, eps, ipg

################################################################################
################################################################################
################################################################################
################################################################################

# compute the residual error 
# where to look for the files?

directoriesCompute = [ \
               'case02_p01_q01_nelem002048_steps00020' \
                ]
                
directoriesTex = ['case02_p01_q01_nelem000032_steps00002', \
               'case02_p01_q01_nelem000128_steps00005', \
               'case02_p01_q01_nelem000512_steps00010', \
               'case02_p01_q01_nelem002048_steps00020' \
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
   
   







