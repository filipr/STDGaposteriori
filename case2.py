from dolfin import *
import numpy as np 

from oneStep import readRHS, dualError
from tex import makeTex
from resError import compute 

directoriesCompute = [ \
               'case02_p02_q01_nelem002048_steps00032' \
               ]
                
directoriesTex = [ \
               'case02_p01_q01_nelem000128_steps00005', \
               'case02_p01_q01_nelem000512_steps00008', \
               'case02_p01_q01_nelem002048_steps00010', \
               'case02_p02_q01_nelem000128_steps00008', \
               'case02_p02_q01_nelem000512_steps00016', \
               'case02_p02_q01_nelem002048_steps00032' \
                ]
#                'case02_p01_q01_nelem008192_steps00020', \


# 'case02_p03_q02_nelem000128_steps00008', \
# 'case02_p03_q02_nelem000512_steps00016', \
# 'case02_p03_q02_nelem002048_steps00032', \

# Print EOC in the tables? 
orders = True

# Compute the dual error OR set the tables with zero value 
computeDual = True
  
   
#Compute the errors and save them to files 
compute(directoriesCompute, computeDual) 

# Make tex table from the given files !!!
makeTex(directoriesTex, orders )
   
   







