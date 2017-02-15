#from dolfin import *
#import numpy as np 

from tex import makeTex
from resError import compute 


directoriesCompute = [ \
               'case24_p02_q01_nelem000288_steps00020' \
                ]
                
directoriesTex = [ \
                 'case24_p02_q01_nelem000288_steps00020' \
                 ]

                ]
# Print EOC in the tables? 
orders = True
  

#Compute the errors and save them to files 
compute(directoriesCompute) 

# Make tex table from the given files !!!
makeTex(directoriesTex, orders)
   
   







