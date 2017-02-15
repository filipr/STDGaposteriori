import numpy as np 


def setTau( p, q, m, h0, tau0): 
   #h0 =  pow(tau0,(q+1.0)/p)
   h = [h0] 
   #tau0 = pow(h0,p/(q+1.0))
   tau = [ tau0 ] 
   
   for i in range(1,m): 
      h.append( h0*pow(2.0,-i) )
      tau.append( tau0*pow(2.0, (-i*p)/(q+1.0) ) )
   
   return h,tau 
   
print 'h & tau = ' , setTau( 2, 1, 3, np.sqrt(2)/4.0, 0.125 ) 

