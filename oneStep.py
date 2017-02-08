from dolfin import *
import numpy as np
import time

def oneStep( fileName ): 
   '''
   compute the residuum for one time step of the time-dependent convection diffusion problem:
   solves the (weighted) poisson problem: 
   -div ( A \grad u ) = g, where g is the residuum given by: 
   g(v) = (f - u',v) - (\grad u, \grad v) - {u}_{m-1} v_{m-1}
   on the domain (0,1)x(0,1)x(t_{m-1}, t_m)
   and boundary conditions: 
   
   the fileName file is exported fro ADGFEM and it should contain: 
   information about the geometry and the RHS for the problem 
   xCoord = number of points with respect to space variables 
   coarseH = diam of elems in the original mesh in adgfem 
   tCoord = number of points with respect to time variables  
   start, ending = t_m-1, t_m 
   dK = parameter of the estimates in adgfem   
   '''

   # read the data
   with open(fileName) as lines: 
       lbCorner_x, lbCorner_y, rtCorner_x, rtCorner_y = \
       [float(x) for x in next(lines).split()] # read the first line 
       xCoord, coarseH, tCoord, start, ending, dK = \
       [float(x) for x in next(lines).split()] # read the 2nd line 
       data = []
       for line in lines:
         data.append([float(v) for v in line.split()]) 
       lines.close()      
      
   # to numpy array : 
   #vals = np.array( [data[i][3] for i in range(len(data))] )  
   tau = ending - start  
   
   coords = np.array( [data[i][0:3] for i in range(len(data))] ) 
   #print 'tau =' , tau
   #print 'h = ' , coarseH 
   #print 'dK = ', dK

   fMinusDer = np.array( [data[i][3] for i in range(len(data))] )
   gradX = np.array( [data[i][4] for i in range(len(data))] )
   gradY = np.array( [data[i][5] for i in range(len(data))] )
   jump = np.array( [data[i][6] for i in range(len(data))] )

   x0 = lbCorner_x
   y0 = lbCorner_y
   z0 = start 
   x1 = rtCorner_x
   y1 = rtCorner_y
   z1 = ending 
   lbCorner = Point( x0, y0, z0 ) 
   rtCorner = Point( x1, y1, z1 )  
   #print 'dsad:' , lbCorner_x, lbCorner_y, rtCorner_x, rtCorner_y
   #print 'startEnd', start, ending

   xCoord = int(xCoord)
   tCoord = int(tCoord)
   nx = xCoord
   ny = xCoord
   nz = tCoord

   Ax = pow(coarseH / dK, 2.) #1.0 #(x1 - x0) / nx
   Ay = pow(coarseH / dK, 2.) #1.0 # (x1 - x0) / ny 
   Az = pow(tau / dK, 2.) #1.0 #(x1 - x0) / nz

   mesh = BoxMesh(lbCorner, rtCorner, nx-1, ny-1, nz-1)
   #print( mesh.coordinates()  )  
   

   d = mesh.geometry().dim()
   print( 'Nelem:',  mesh.num_entities(0), mesh.num_entities(3)  ) 
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
                   ('0.0','0.0', 'C')), A = Ax, B = Ay, C = Az)
   C = as_matrix(c)

   # Define variational problem
   u = TrialFunction(V)
   v = TestFunction(V)

   a = inner(C*grad(u), grad(v))*dx
   L = f*v*dx + dWdx*v.dx(0)*dx + dWdy*v.dx(1)*dx + wJump*v*dsJump

   # Compute solution
#   solver = KrylovSolver("gmres", "ilu") # amg, ilu, schwarz
#   solver.parameters["relative_tolerance"] = tol
#   solver.parameters["absolute_tolerance"] = tol
#   solver.parameters["maximum_iterations"] = 4000
#   solver.parameters['gmres']['restart'] = 40
#   # do not start with zero vector,but some initial guess
#   solver.parameters['nonzero_initial_guess'] = False
#   solver.parameters["monitor_convergence"] = False  
#   A = assemble(a) 
#   F = assemble(L)
#   #print 'F=', type(F), F.array()
#   krylovIter = solver.solve( A, u.vector(), F )  
   u = Function(V)
   problem = LinearVariationalProblem(a, L, u, bc)
   solver = LinearVariationalSolver(problem) 
   solver.parameters["linear_solver"] = "gmres"
   solver.parameters["preconditioner"] = "ilu"
   #solver.parameters["krylov_solver"]["monitor_convergence"] = True
   start = time.clock()
   solver.solve()
   print 'Time elapsed:' , (time.clock() - start)
   print 'Zbytecne se sestavuje system znovu? '
   
   #solve(a == L, u, bc,  \
   #         solver_parameters={'linear_solver': 'cg', \
   #                      'preconditioner': 'ilu', 'monitor_convergence': True}) # icc, ilu, amg, sor   

   # estimate the error
   zero = Constant( 0.0 ) 
   zeroFun = interpolate(zero, V)


#   errorL2 = errornorm( u, zeroFun, norm_type='L2', degree_rise=1)
#   print 'E L2: ' , errorL2
#   ## WATCH it computes H1 NORM not seminorm !
#   errorH1 = errornorm( u, zeroFun, norm_type='H1') 
#   print 'E_H1 = ' , errorH1  
   
   # weighted H1 seminorm SQUARED
   M = ( Ax*u.dx(0)**2. + Ay*u.dx(1)**2. + Az*u.dx(2)**2. ) *dx
   error = assemble(M) 
   
   return error 
