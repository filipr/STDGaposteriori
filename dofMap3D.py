from dolfin import *
import numpy as np

# we have the values of a functio in the mesh coor

# read the input file with the coordinates and function values

# test in3D.txt 
i = 1
fileName = '/home/filip/adgfem/RTN/rtn_0' + str(i)  + '.txt' 
print(fileName) 


with open(fileName) as lines: 
    xCoord, coarseH, tCoord, start, ending, dK = [float(x) for x in next(lines).split()] # read first line
    data = []
    for line in lines:
      data.append([float(v) for v in line.split()]) 
      
# to numpy array : 
#vals = np.array( [data[i][3] for i in range(len(data))] )  
tau = ending - start 

coords = np.array( [data[i][0:3] for i in range(len(data))] ) 



fMinusDer = np.array( [data[i][3] for i in range(len(data))] )
gradX = np.array( [data[i][4] for i in range(len(data))] )
gradY = np.array( [data[i][5] for i in range(len(data))] )
jump = np.array( [data[i][6] for i in range(len(data))] )


#print( 'jump=', jump) 


#mesh = UnitCubeMesh(2, 2, 2) 
x0 = 0.0 
y0 = 0.0 
z0 = start 
x1 = 1.0 
y1 = 1.0 
z1 = ending 
lbCorner = Point( x0, y0, z0 ) 
rtCorner = Point( x1, y1, z1 ) 

xCoord = int(xCoord)
tCoord = int(tCoord)
nx = xCoord
ny = xCoord
nz = tCoord

hx = coarseH #1.0 #(x1 - x0) / nx
hy = coarseH #1.0 # (x1 - x0) / ny 
hz = tau #1.0 #(x1 - x0) / nz

mesh = BoxMesh(lbCorner, rtCorner, nx-1, ny-1, nz-1)
#print( mesh.coordinates()  )  

d = mesh.geometry().dim()
#center = Point( 0.5, 0.5, 0.7 )

V = FunctionSpace(mesh, 'CG', 1) 

# Coordinates of all vertices
vertex_x = mesh.coordinates().reshape((-1, d))

# vi_vx[vertex_index] = coordinates of index
vi_vx = vertex_x.tolist() 


if ( ((abs(vi_vx - coords)) > 1.E-9).any() ):
   print( 'Wrong coordinates!' ) 
   quit() 
   
   
# f - u'
f = Function(V) 
dWdx = Function(V) 
dWdy = Function(V) 
wJump = Function(V) 
#v_d = vertex_to_dof_map(V)
d_v = dof_to_vertex_map(V)

for i in range( len(fMinusDer) ): 
   #print( i, v_d[i] )
   f.vector()[i] = fMinusDer[d_v[i]] 
   dWdx.vector()[i] = gradX[d_v[i]] 
   dWdy.vector()[i] = gradY[d_v[i]] 
   wJump.vector()[i] = jump[d_v[i]] 
   

# Define Dirichlet boundary (x = 0 or x = 1) or y=0 or y=1.0
def boundary(x):
    return x[0] < x0 + DOLFIN_EPS or x[0] > x1 - DOLFIN_EPS or \
           x[1] < y0 + DOLFIN_EPS or x[1] > y1 - DOLFIN_EPS

# mark bottom for the jump
boundary_parts = FacetFunction('size_t', mesh)  
bottom = AutoSubDomain(lambda x: near(x[2], start )) 
bottom.mark(boundary_parts, 1)
dsJump = Measure("ds", subdomain_id = 1, subdomain_data = boundary_parts)


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
L = f*v*dx + dWdx*v.dx(0)*dx + dWdy*v.dx(1)*dx + wJump*v*dsJump


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





