#from dolfin import *
#import numpy as np 

from tex import makeTex
from resError import compute 


directoriesCompute = [ \
               'case24_p02_q01_nelem004608_steps00040' \
                ]
                
directoriesTex = [ \
               'case24_p02_q01_nelem000288_steps00010', \
               'case24_p02_q01_nelem001152_steps00020', \
               'case24_p02_q01_nelem004608_steps00040' \
                 ]

# case24_p02_q01_nelem004608_steps00040


# Print EOC in the tables? 
orders = True

# Compute the dual error OR set the tables with zero value 
computeDual = False
  
   
#Compute the errors and save them to files 
compute(directoriesCompute, computeDual) 

# Make tex table from the given files !!!
makeTex(directoriesTex, orders )





