from dolfin import *
import numpy as np


mesh = UnitCubeMesh(2, 2, 2)
d = mesh.geometry().dim()


V = FunctionSpace(mesh, 'CG', 1)

# Coordinates of all vertices
vertex_x = mesh.coordinates().reshape((-1, d))

# vi_vx[vertex_index] = coordinates of index
vi_vx = vertex_x.tolist()
#print( 'vi_dx:', vi_vx )
 
with open( 'in3D.txt', 'w') as f:
   for v in vi_vx: 
      f.write( '{:2} {:} {:2} {:2} {:1}'.format( v[0], v[1], v[2], 1.0,  '\n' ) ) 
   f.close()


