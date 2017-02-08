from dolfin import *
import numpy as np

# we have the values of a functio in the mesh coor

# read the input file with the coordinates and function values
data = []
with open('in.txt') as lines:
   for line in lines:
      data.append([float(v) for v in line.split()]) 
      
# to numpy array : 
vals = np.array(  [data[i][2] for i in range(len(data))] ) 


print( 'x=', vals ) 


#for i in range(len(data)):
#   print( 'x[',i,'] = ', data[i] )


mesh = UnitSquareMesh(2, 2)
d = mesh.geometry().dim()


V = FunctionSpace(mesh, 'CG', 1)


f = Function(V) 
v_d = vertex_to_dof_map(V)

d_v = dof_to_vertex_map(V)

for i in range( len(vals) ): 
   print( i, v_d[i] )
   f.vector()[i] = vals[d_v[i]]
plot(f, interactive = True) 



v_d = vertex_to_dof_map(V)
print( 'vd:', v_d )
d_v = dof_to_vertex_map(V)
print( 'd_v:', d_v )
# Coordinates of all dofs
dof_x = V.dofmap().tabulate_all_coordinates(mesh).reshape((-1, d))

# di_dx[dof_index] = coordinates of dof
di_dx = dof_x.tolist()
print( 'di_dx:', di_dx )
## Coordinates of all vertices
#vertex_x = mesh.coordinates().reshape((-1, d))

## vi_vx[vertex_index] = coordinates of index
#vi_vx = vertex_x.tolist()
#print( 'vi_dx:', vi_vx )

#n = V.dofmap().num_entity_dofs(0)
## The test shows how the vertex_to_dofmap and dof_to_vertex_map
## should be interpreted. [vi/n] reflects the shift due to mixed space

#assert all(np.allclose(di_dx[di], vi_vx[vi/n])
#               for vi, di in enumerate(v_d))
#assert all(np.allclose(di_dx[di], vi_vx[int(vi)/n])
#               for di, vi in enumerate(d_v))
