from dolfin import *
import numpy as np

# we have the values of a functio in the mesh coor

# read the input file with the coordinates and function values
data = []
with open('in3D.txt') as lines:
   for line in lines:
      data.append([float(v) for v in line.split()]) 
      
# to numpy array : 
vals = np.array(  [data[i][3] for i in range(len(data))] ) 

#print( 'x=', vals ) 

#mesh = UnitCubeMesh(2, 2, 2) 
x0 = 0.0 
y0 = 0.0 
z0 = 0.0 
x1 = 1.0 
y1 = 1.0 
z1 = 1.0 
lbCorner = Point( x0, y0, z0 ) 
rtCorner = Point( x1, y1, z1 ) 

nx = 2 
ny = 2 
nz = 2 

hx = (x1 - x0) / nx
hy = (x1 - x0) / ny 
hz = (x1 - x0) / nz

mesh = BoxMesh(lbCorner, rtCorner, nx, ny, nz)
d = mesh.geometry().dim()
center = Point( 0.5, 0.5, 0.7 )

V = FunctionSpace(mesh, 'CG', 1)


f = Function(V) 
v_d = vertex_to_dof_map(V)
d_v = dof_to_vertex_map(V)

for i in range( len(vals) ): 
   #print( i, v_d[i] )
   f.vector()[i] = vals[d_v[i]]
   
print( 'f(center):', f(center) ) 


# Define Dirichlet boundary (x = 0 or x = 1) or y=0 or y=1.0
def boundary(x):
    return x[0] < x0 + DOLFIN_EPS or x[0] > x1 - DOLFIN_EPS or \
           x[1] < y0 + DOLFIN_EPS or x[1] > y1 - DOLFIN_EPS
# Define boundary condition
u0 = Constant(0.0)
bc = DirichletBC(V, u0, boundary) 

c = Expression((('A','0.0', '0.0'),
                ('0.0', 'B', '0.0'),
                ('0.0','0.0', 'C')), A = hx, B = hy, C = hz)
C = as_matrix(c)

# Define variational problem
u = TrialFunction(V)
v = TestFunction(V)

a = inner(C*grad(u), grad(v))*dx
L = f*v*dx

# Compute solution
u = Function(V)
solve(a == L, u, bc) 

zero = Constant( 0.0 ) 
zeroFun = interpolate(zero, V)


errorL2 = errornorm( u, zeroFun, norm_type='L2', degree_rise=1)
print 'E L2: ' , errorL2
## WATCH it computes H1 NORM not seminorm !
errorH1 = errornorm( u, zeroFun, norm_type='H1') 
print 'E_H1 = ' , errorH1 





