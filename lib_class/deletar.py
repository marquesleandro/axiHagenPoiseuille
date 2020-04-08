
class AxiElement2D:
 def __init__(_self, _z, _r, _IEN, _GAUSSPOINTS):
  _self.z = _z
  _self.r = _r
  _self.IEN = _IEN


  if _GAUSSPOINTS == 3:
   _self.NUMGAUSS = 3  #Number of Gauss Points


   #                                 l1                 l2
   _self.GQPoints = np.array([[0.16666666666667, 0.16666666666667], 
                              [0.16666666666667, 0.66666666666667], 
                              [0.66666666666667, 0.16666666666667]])


   #                                    w
   _self.GQWeights = np.array([[0.333333333333333], 
                               [0.333333333333333], 
                               [0.333333333333333]])





  elif _GAUSSPOINTS == 4:
   _self.NUMGAUSS = 4  #Number of Gauss Points


   #                                 l1                 l2
   _self.GQPoints = np.array([[0.33333333333333, 0.33333333333333], 
                              [0.20000000000000, 0.20000000000000], 
                              [0.20000000000000, 0.60000000000000], 
                              [0.60000000000000, 0.20000000000000]])


   #                                    w
   _self.GQWeights = np.array([[-0.56250000000000], 
                               [0.520833333333333], 
                               [0.520833333333333], 
                               [0.520833333333333]])



  
  elif _GAUSSPOINTS == 6:
   _self.NUMGAUSS = 6  #Number of Gauss Points


   #                                 l1                 l2
   _self.GQPoints = np.array([[0.44594849091597, 0.44594849091597],
                              [0.44594849091597, 0.10810301816807],
                              [0.10810301816807, 0.44594849091597],
                              [0.09157621350977, 0.09157621350977],
                              [0.09157621350977, 0.81684757298046],
                              [0.81684757298046, 0.09157621350977]])


   #                                    w
   _self.GQWeights = np.array([[0.22338158967801],
                               [0.22338158967801],
                               [0.22338158967801],
                               [0.10995174365532],
                               [0.10995174365532],
                               [0.10995174365532]])


  elif _GAUSSPOINTS == 12:
   _self.NUMGAUSS = 12  #Number of Gauss Points

   #                                 l1                 l2
   _self.GQPoints = np.array([[0.24928674517091, 0.24928674517091],
                              [0.24928674517091, 0.50142650965818],
                              [0.50142650965818, 0.24928674517091],
                              [0.06308901449150, 0.06308901449150],
                              [0.06308901449150, 0.87382197101700],
                              [0.87382197101700, 0.06308901449150],
                              [0.31035245103378, 0.63650249912140],
                              [0.63650249912140, 0.05314504984482],
                              [0.05314504984482, 0.31035245103378],
                              [0.63650249912140, 0.31035245103378],
                              [0.31035245103378, 0.05314504984482],
                              [0.05314504984482, 0.63650249912140]])

   #                                    w
   _self.GQWeights = np.array([[0.11678627572638],
                               [0.11678627572638],
                               [0.11678627572638],
                               [0.05084490637021],
                               [0.05084490637021],
                               [0.05084490637021],
                               [0.08285107561837],
                               [0.08285107561837],
                               [0.08285107561837],
                               [0.08285107561837],
                               [0.08285107561837],
                               [0.08285107561837]])
   '''
 
 
    #                                 l1                 l2
   _self.GQPoints = np.array([[0.219429982550000, 0.561140034900000],
                              [0.561140034900000, 0.219429982550000],
                              [0.219429982550000, 0.219429982550000],
                              [0.480137964112000, 0.039724071775600],
                              [0.039724071775600, 0.480137964112000],
                              [0.480137964112000, 0.480137964112000],
                              [0.019371724361200, 0.839009259715000],
                              [0.141619015924000, 0.019371724361200],
                              [0.839009259715000, 0.141619015924000],
                              [0.141619015924000, 0.839009259715000],
                              [0.019371724361200, 0.141619015924000],
                              [0.839009259715000, 0.019371724361200]])


    #                                    w
   _self.GQWeights = np.array([[0.171333124153000],
                               [0.171333124153000],
                               [0.171333124153000],
                               [0.080731089593000],
                               [0.080731089593000],
                               [0.080731089593000],
                               [0.040634559793700],
                               [0.040634559793700],
                               [0.040634559793700],
                               [0.040634559793700],
                               [0.040634559793700],
                               [0.040634559793700]])
 
 
  '''
  else:
   print ""
   print " Error: Gauss Points not found"
   print ""
   sys.exit()


 def linear(_self,_e):
  _self.NUMNODE = 3  #Linear Triangular Element - 3 Nodes
  
  N = np.zeros([_self.NUMGAUSS,_self.NUMNODE], dtype = float)

  dNdl1 = np.zeros([_self.NUMGAUSS,_self.NUMNODE], dtype = float)
  dNdl2 = np.zeros([_self.NUMGAUSS,_self.NUMNODE], dtype = float)
  
  dzdl1 = np.zeros([_self.NUMGAUSS,1], dtype = float)
  dzdl2 = np.zeros([_self.NUMGAUSS,1], dtype = float)
  drdl1 = np.zeros([_self.NUMGAUSS,1], dtype = float)
  drdl2 = np.zeros([_self.NUMGAUSS,1], dtype = float)

  dNdz = np.zeros([_self.NUMGAUSS,_self.NUMNODE], dtype = float)
  dNdr = np.zeros([_self.NUMGAUSS,_self.NUMNODE], dtype = float)

  J = np.zeros([2,2], dtype = float)
  jacobian = np.zeros([_self.NUMGAUSS,1], dtype = float)


  # A loop is required for each pair of coordinates (l1,l2)
  for k in range(0,_self.NUMGAUSS):
    
   # Area Coordinates
   # Lewis pag. 67 Eq. 3.129
   L1 = 1.0 - _self.GQPoints[k][0] - _self.GQPoints[k][1]   #L1 = 1 - l1 - l2
   L2 = _self.GQPoints[k][0]                                #L2 = l1
   L3 = _self.GQPoints[k][1]                                #L3 = l2

   # Shape Functions
   # Lewis pag. 67 Eq. 3.129
   N[k][0] = L1  #N1 = L1
   N[k][1] = L2  #N2 = L2
   N[k][2] = L3  #N3 = L3

   # Shape Functions Derivatives in respect to l1
   dNdl1[k][0] = -1.0   #dN1/dl1
   dNdl1[k][1] =  1.0   #dN2/dl1
   dNdl1[k][2] =  0.0   #dN3/dl1

   # Shape Functions Derivatives in respect to l2
   dNdl2[k][0] = -1.0   #dN1/dl2
   dNdl2[k][1] =  0.0   #dN2/dl2
   dNdl2[k][2] =  1.0   #dN3/dl2

   # Coordinate Transfomation
   # Lewis pag. 64 Eq. 3.108 for 1D
   # Lewis pag. 66 Eq. 3.121 for 2D
   for i in range(0,_self.NUMNODE):
    ii = _self.IEN[_e][i]
    
    dzdl1[k] += _self.z[ii]*dNdl1[k][i]   # dz/dl1 = z1*dN1/dl1 + z2*dN2/dl1 + z3*dN3/dl1 ...
    drdl1[k] += _self.r[ii]*dNdl1[k][i]   # dr/dl1 = r1*dN1/dl1 + r2*dN2/dl1 + r3*dN3/dl1 ...
    dzdl2[k] += _self.z[ii]*dNdl2[k][i]   # dz/dl2 = z1*dN1/dl2 + z2*dN2/dl2 + z3*dN3/dl2 ...
    drdl2[k] += _self.r[ii]*dNdl2[k][i]   # dr/dl2 = r1*dN1/dl2 + r2*dN2/dl2 + r3*dN3/dl2 ...


   # Jacobian Matrix
   # Lewis pag. 66 Eq. 3.121
   J[0][0] = dzdl1[k]
   J[0][1] = drdl1[k]
   J[1][0] = dzdl2[k]
   J[1][1] = drdl2[k]

   jacobian[k] = np.linalg.det(J)


   # Shape Functions Derivatives in respect to x and y
   # Lewis pag. 65 Eq. 3.116
   for i in range(0,_self.NUMNODE):
    dNdz[k][i] = (1.0/jacobian[k])*( dNdl1[k][i]*drdl2[k] - dNdl2[k][i]*drdl1[k])
    dNdr[k][i] = (1.0/jacobian[k])*(-dNdl1[k][i]*dzdl2[k] + dNdl2[k][i]*dzdl1[k])



  _self.mass = np.zeros([_self.NUMNODE,_self.NUMNODE], dtype = float)
  _self.kzz = np.zeros([_self.NUMNODE,_self.NUMNODE], dtype = float)
  _self.kzr = np.zeros([_self.NUMNODE,_self.NUMNODE], dtype = float)
  _self.krz = np.zeros([_self.NUMNODE,_self.NUMNODE], dtype = float)
  _self.krr = np.zeros([_self.NUMNODE,_self.NUMNODE], dtype = float)
  _self.gz = np.zeros([_self.NUMNODE,_self.NUMNODE], dtype = float)
  _self.gr = np.zeros([_self.NUMNODE,_self.NUMNODE], dtype = float)
  _self.dz = np.zeros([_self.NUMNODE,_self.NUMNODE], dtype = float)
  _self.dr = np.zeros([_self.NUMNODE,_self.NUMNODE], dtype = float)
  

  # Elementary Matrices
  # P.S: It is divided by 2 due to relation 
  # of parallelogram and triangle areas --> DxDy = (jacobian*Dl1*Dl2)/2
  for k in range(0,_self.NUMGAUSS): 
   for i in range(0,_self.NUMNODE):
    for j in range(0,_self.NUMNODE):
    
     _self.mass[i][j] += N[k][i]*N[k][j]*jacobian[k]*_self.GQWeights[k]/2.0
    
     _self.kzz[i][j] += dNdz[k][i]*dNdz[k][j]*jacobian[k]*_self.GQWeights[k]/2.0
     _self.kzr[i][j] += dNdz[k][i]*dNdr[k][j]*jacobian[k]*_self.GQWeights[k]/2.0
     _self.krz[i][j] += dNdr[k][i]*dNdz[k][j]*jacobian[k]*_self.GQWeights[k]/2.0
     _self.krr[i][j] += dNdr[k][i]*dNdr[k][j]*jacobian[k]*_self.GQWeights[k]/2.0
    
     _self.gz[i][j] += dNdz[k][j]*N[k][i]*jacobian[k]*_self.GQWeights[k]/2.0
     _self.gr[i][j] += dNdr[k][j]*N[k][i]*jacobian[k]*_self.GQWeights[k]/2.0
    
     _self.dz[i][j] += dNdz[k][j]*N[k][i]*jacobian[k]*_self.GQWeights[k]/2.0
     _self.dr[i][j] += dNdr[k][j]*N[k][i]*jacobian[k]*_self.GQWeights[k]/2.0
  

 def mini(_self,_e):
  _self.NUMNODE = 4  #Mini Triangular Element - 4 Nodes
  
  N = np.zeros([_self.NUMGAUSS,_self.NUMNODE], dtype = float)

  dNdl1 = np.zeros([_self.NUMGAUSS,_self.NUMNODE], dtype = float)
  dNdl2 = np.zeros([_self.NUMGAUSS,_self.NUMNODE], dtype = float)
  
  dzdl1 = np.zeros([_self.NUMGAUSS,1], dtype = float)
  dzdl2 = np.zeros([_self.NUMGAUSS,1], dtype = float)
  drdl1 = np.zeros([_self.NUMGAUSS,1], dtype = float)
  drdl2 = np.zeros([_self.NUMGAUSS,1], dtype = float)

  dNdz = np.zeros([_self.NUMGAUSS,_self.NUMNODE], dtype = float)
  dNdr = np.zeros([_self.NUMGAUSS,_self.NUMNODE], dtype = float)

  J = np.zeros([2,2], dtype = float)
  jacobian = np.zeros([_self.NUMGAUSS,1], dtype = float)



  # A loop is required for each pair of coordinates (l1,l2)
  for k in range(0,_self.NUMGAUSS):
    
   # Area Coordinates
   L1 = _self.GQPoints[k][0]                                #L1 = l1
   L2 = _self.GQPoints[k][1]                                #L2 = l2
   L3 = 1.0 - _self.GQPoints[k][0] - _self.GQPoints[k][1]   #L3 = 1 - l1 - l2
   
   # Shape Functions
   # Mini2D matlab Gustavo
   N[k][0] = L1 - 9.0*L1*L2*L3    #N1 = L1-9*L1*L2*L3
   N[k][1] = L2 - 9.0*L1*L2*L3    #N2 = L2-9*L1*L2*L3
   N[k][2] = L3 - 9.0*L1*L2*L3    #N3 = L3-9*L1*L2*L3
   N[k][3] = 27.0*L1*L2*L3        #N4 = 27*L1*L2*L3

   # Shape Functions Derivatives in respect to l1
   #dN1/dl1
   dNdl1[k][0] =  18.0*_self.GQPoints[k][0]*_self.GQPoints[k][1] + 9.0*_self.GQPoints[k][1]**2\
                 - 9.0*_self.GQPoints[k][1] + 1.0 

   #dN2/dl1
   dNdl1[k][1] =  18.0*_self.GQPoints[k][0]*_self.GQPoints[k][1] + 9.0*_self.GQPoints[k][1]**2\
                 - 9.0*_self.GQPoints[k][1] 

   #dN3/dl1
   dNdl1[k][2] =  18.0*_self.GQPoints[k][0]*_self.GQPoints[k][1] + 9.0*_self.GQPoints[k][1]**2\
                 - 9.0*_self.GQPoints[k][1] - 1.0

   #dN4/dl1
   dNdl1[k][3] =  27.0*(- 2.0*_self.GQPoints[k][0]*_self.GQPoints[k][1] - _self.GQPoints[k][1]**2\
                 + _self.GQPoints[k][1])


   # Shape Functions Derivatives in respect to l2
   #dN1/dl2
   dNdl2[k][0] =  18.0*_self.GQPoints[k][0]*_self.GQPoints[k][1] + 9.0*_self.GQPoints[k][0]**2\
                 - 9.0*_self.GQPoints[k][0] 

   #dN2/dl2
   dNdl2[k][1] =  18.0*_self.GQPoints[k][0]*_self.GQPoints[k][1] + 9.0*_self.GQPoints[k][0]**2\
                 - 9.0*_self.GQPoints[k][0] + 1.0 

   #dN3/dl2
   dNdl2[k][2] =  18.0*_self.GQPoints[k][0]*_self.GQPoints[k][1] + 9.0*_self.GQPoints[k][0]**2\
                 - 9.0*_self.GQPoints[k][0] - 1.0

   #dN4/dl2
   dNdl2[k][3] =  27.0*(- 2.0*_self.GQPoints[k][0]*_self.GQPoints[k][1] - _self.GQPoints[k][0]**2\
                 + _self.GQPoints[k][0])


   # Coordinate Transfomation
   # Lewis pag. 64 Eq. 3.108 for 1D
   # Lewis pag. 66 Eq. 3.121 for 2D
   for i in range(0,_self.NUMNODE):
    ii = _self.IEN[_e][i]
    
    dzdl1[k] += _self.z[ii]*dNdl1[k][i]   # dx/dl1 = x1*dN1/dl1 + x2*dN2/dl1 + x3*dN3/dl1 ...
    drdl1[k] += _self.r[ii]*dNdl1[k][i]   # dy/dl1 = y1*dN1/dl1 + y2*dN2/dl1 + y3*dN3/dl1 ...
    dzdl2[k] += _self.z[ii]*dNdl2[k][i]   # dx/dl2 = x1*dN1/dl2 + x2*dN2/dl2 + x3*dN3/dl2 ...
    drdl2[k] += _self.r[ii]*dNdl2[k][i]   # dy/dl2 = y1*dN1/dl2 + y2*dN2/dl2 + y3*dN3/dl2 ...

   # Jacobian Matrix
   # Lewis pag. 66 Eq. 3.121
   J[0][0] = dzdl1[k]
   J[0][1] = drdl1[k]
   J[1][0] = dzdl2[k]
   J[1][1] = drdl2[k]

   jacobian[k] = np.linalg.det(J)


   # Shape Functions Derivatives in respect to x and y
   # Lewis pag. 65 Eq. 3.116
   for i in range(0,_self.NUMNODE):
    dNdz[k][i] = (1.0/jacobian[k])*( dNdl1[k][i]*drdl2[k] - dNdl2[k][i]*drdl1[k])
    dNdr[k][i] = (1.0/jacobian[k])*(-dNdl1[k][i]*dzdl2[k] + dNdl2[k][i]*dzdl1[k])



  _self.mass = np.zeros([_self.NUMNODE,_self.NUMNODE], dtype = float)
  _self.kzz = np.zeros([_self.NUMNODE,_self.NUMNODE], dtype = float)
  _self.kzr = np.zeros([_self.NUMNODE,_self.NUMNODE], dtype = float)
  _self.krz = np.zeros([_self.NUMNODE,_self.NUMNODE], dtype = float)
  _self.krr = np.zeros([_self.NUMNODE,_self.NUMNODE], dtype = float)
  _self.gz = np.zeros([_self.NUMNODE,_self.NUMNODE], dtype = float)
  _self.gr = np.zeros([_self.NUMNODE,_self.NUMNODE], dtype = float)
  _self.dz = np.zeros([_self.NUMNODE,_self.NUMNODE], dtype = float)
  _self.dr = np.zeros([_self.NUMNODE,_self.NUMNODE], dtype = float)
  

  # Elementary Matrices
  # P.S: It is divided by 2 due to relation 
  # of parallelogram and triangle areas --> DxDy = (jacobian*Dl1*Dl2)/2
  for k in range(0,_self.NUMGAUSS): 
   for i in range(0,_self.NUMNODE):
    for j in range(0,_self.NUMNODE):
    
     _self.mass[i][j] += N[k][i]*N[k][j]*jacobian[k]*_self.GQWeights[k]/2.0
    
     _self.kzz[i][j] += dNdz[k][i]*dNdz[k][j]*jacobian[k]*_self.GQWeights[k]/2.0
     _self.kzr[i][j] += dNdz[k][i]*dNdr[k][j]*jacobian[k]*_self.GQWeights[k]/2.0
     _self.krz[i][j] += dNdr[k][i]*dNdz[k][j]*jacobian[k]*_self.GQWeights[k]/2.0
     _self.krr[i][j] += dNdr[k][i]*dNdr[k][j]*jacobian[k]*_self.GQWeights[k]/2.0
    
     _self.gz[i][j] += dNdz[k][j]*N[k][i]*jacobian[k]*_self.GQWeights[k]/2.0
     _self.gr[i][j] += dNdr[k][j]*N[k][i]*jacobian[k]*_self.GQWeights[k]/2.0
    
     _self.dz[i][j] += dNdz[k][j]*N[k][i]*jacobian[k]*_self.GQWeights[k]/2.0
     _self.dr[i][j] += dNdr[k][j]*N[k][i]*jacobian[k]*_self.GQWeights[k]/2.0
  





 def quadratic(_self,_e):
  _self.NUMNODE = 6  #Quadratic Triangular Element - 6 Nodes
  
#  _self.N = np.zeros([_self.NUMGAUSS,_self.NUMNODE], dtype = float)

#  _self.dNdl1 = np.zeros([_self.NUMGAUSS,_self.NUMNODE], dtype = float)
#  _self.dNdl2 = np.zeros([_self.NUMGAUSS,_self.NUMNODE], dtype = float)
  
  dzdl1 = np.zeros([_self.NUMGAUSS,1], dtype = float)
  dzdl2 = np.zeros([_self.NUMGAUSS,1], dtype = float)
  drdl1 = np.zeros([_self.NUMGAUSS,1], dtype = float)
  drdl2 = np.zeros([_self.NUMGAUSS,1], dtype = float)

  dNdz = np.zeros([_self.NUMGAUSS,_self.NUMNODE], dtype = float)
  dNdr = np.zeros([_self.NUMGAUSS,_self.NUMNODE], dtype = float)

  J = np.zeros([2,2], dtype = float)
  jacobian = np.zeros([_self.NUMGAUSS,1], dtype = float)


  _self.N = np.array([[-0.124998982535, -0.124998982535, 0.001430579518, 0.248575525272, 0.499995930140, 0.499995930140],
                      [-0.124998982535, 0.001430579518, -0.124998982535, 0.499995930140, 0.499995930140, 0.248575525272],
                      [0.001430579518, -0.124998982535, -0.124998982535, 0.499995930140, 0.248575525272, 0.499995930140],
                      [-0.055128566992, -0.055128566992, 0.653307703047, 0.015920894998, 0.220514267970, 0.220514267970],
                      [-0.055128566992, 0.653307703047, -0.055128566992, 0.220514267970, 0.220514267970, 0.015920894998],
                      [0.653307703047, -0.055128566992, -0.055128566992, 0.220514267970, 0.015920894998, 0.220514267970],
                      [-0.117715163308, 0.173768363654, -0.047496257199, 0.790160442766, 0.135307828169, 0.065974785919],
                      [0.173768363654, -0.047496257199, -0.117715163308, 0.135307828169, 0.065974785919, 0.790160442766],
                      [-0.047496257199, -0.117715163308, 0.173768363654, 0.065974785919, 0.790160442766, 0.135307828169],
                      [0.173768363654, -0.117715163308, -0.047496257199, 0.790160442766, 0.065974785919, 0.135307828169],
                      [-0.117715163308, -0.047496257199, 0.173768363654, 0.065974785919, 0.135307828169, 0.790160442766],
                      [-0.047496257199, 0.173768363654, -0.117715163308, 0.135307828169, 0.790160442766, 0.065974785919]])


  _self.dNdl1 = np.array([[-0.002853019316, 0.000000000000, -1.005706038632, 0.997146980684, -0.997146980684, 1.008559057948],
                          [-0.002853019316, 0.000000000000, 0.002853019316, 2.005706038632, -2.005706038632, 0.000000000000 ],
                          [1.005706038632, 0.000000000000, 0.002853019316, 0.997146980684,  -0.997146980684, -1.008559057948],
                          [-0.747643942034, 0.000000000000, -2.495287884068, 0.252356057966, -0.252356057966, 3.242931826102],
                          [-0.747643942034, 0.000000000000, 0.747643942034, 3.495287884068,  -3.495287884068, -0.00000000000],
                          [2.495287884068, 0.000000000000, 0.747643942034, 0.252356057966, -0.252356057966, -3.242931826102 ],
                          [0.241409804136, 0.000000000000, 0.787419800620, 2.546009996484, -2.546009996484, -1.028829604756 ],
                          [1.546009996484, 0.000000000000, -0.241409804137, 0.212580199379, -0.212580199379, -1.304600192347],
                          [-0.787419800621, 0.000000000000, -1.546009996485, 1.241409804136, -1.241409804136, 2.333429797106],
                          [1.546009996484, 0.000000000000, 0.787419800620, 1.241409804136, -1.241409804136, -2.333429797104 ],
                          [0.241409804136, 0.000000000000, -1.546009996485, 0.212580199379, -0.212580199379, 1.304600192349 ],
                          [-0.787419800621, 0.000000000000, -0.241409804137, 2.546009996484, -2.546009996484, 1.028829604758]])


  _self.dNdl2 = np.array([[0.000000000000, -0.002853019316, -1.005706038632, 0.997146980684, 1.008559057948, -0.997146980684],
                          [0.000000000000, 1.005706038632, 0.002853019316, 0.997146980684, -1.008559057948, -0.997146980684 ],
                          [0.000000000000, -0.002853019316, 0.002853019316, 2.005706038632, -0.000000000000, -2.005706038632],
                          [0.000000000000, -0.747643942034, -2.495287884068, 0.252356057966, 3.242931826102, -0.252356057966],
                          [0.000000000000, 2.495287884068, 0.747643942034, 0.252356057966, -3.242931826102, -0.252356057966 ],
                          [0.000000000000, -0.747643942034, 0.747643942034, 3.495287884068, -0.000000000000, -3.495287884068],
                          [0.000000000000, 1.546009996484, 0.787419800620, 1.241409804136, -2.333429797104, -1.241409804136 ],
                          [0.000000000000, -0.787419800621, -0.241409804137, 2.546009996484, 1.028829604758, -2.546009996484],
                          [0.000000000000, 0.241409804136, -1.546009996485, 0.212580199379, 1.304600192349, -0.212580199379 ],
                          [0.000000000000, 0.241409804136, 0.787419800620, 2.546009996484, -1.028829604756, -2.546009996484 ],
                          [0.000000000000, -0.787419800621, -1.546009996485, 1.241409804136, 2.333429797106, -1.241409804136],
                          [0.000000000000, 1.546009996484, -0.241409804137, 0.212580199379, -1.304600192347, -0.212580199379]])

  # A loop is required for each pair of coordinates (l1,l2)
  for k in range(0,_self.NUMGAUSS):
    
   # Area Coordinates
   # Lewis pag. 67 Eq. 3.129
#   L1 = 1.0 - _self.GQPoints[k][0] - _self.GQPoints[k][1]   #L1 = 1 - l1 - l2
#   L2 = _self.GQPoints[k][0]                                #L2 = l1
#   L3 = _self.GQPoints[k][1]                                #L3 = l2

   # Shape Functions
   # Lewis pag. 67 Eq. 3.130
#   _self.N[k][0] = L1*(2.0*L1 - 1.0)  #N1 = L1*(2*L1-1)
#   _self.N[k][1] = L2*(2.0*L2 - 1.0)  #N2 = L2*(2*L2-1)
#   _self.N[k][2] = L3*(2.0*L3 - 1.0)  #N3 = L3*(2*L3-1)
#   _self.N[k][3] = 4.0*L1*L2          #N4 = 4*L1*L2
#   _self.N[k][4] = 4.0*L2*L3          #N5 = 4*L2*L3
#   _self.N[k][5] = 4.0*L3*L1          #N6 = 4*L3*L1

   # Shape Functions Derivatives in respect to l1
#   _self.dNdl1[k][0] = -3.0 + 4.0*_self.GQPoints[k][0] + 4.0*_self.GQPoints[k][1]   #dN1/dl1
#   _self.dNdl1[k][1] =  4.0*_self.GQPoints[k][0] - 1.0                              #dN2/dl1
#   _self.dNdl1[k][2] =  0.0                                                         #dN3/dl1
#   _self.dNdl1[k][3] =  4.0 - 8.0*_self.GQPoints[k][0] - 4.0*_self.GQPoints[k][1]   #dN4/dl1
#   _self.dNdl1[k][4] =  4.0*_self.GQPoints[k][1]                                    #dN5/dl1
#   _self.dNdl1[k][5] = -4.0*_self.GQPoints[k][1]                                    #dN6/dl1


   # Shape Functions Derivatives in respect to l2
#   _self.dNdl2[k][0] = -3.0 + 4.0*_self.GQPoints[k][0] + 4.0*_self.GQPoints[k][1]   #dN1/dl2
#   _self.dNdl2[k][1] =  0.0                                                         #dN2/dl2
#   _self.dNdl2[k][2] =  4.0*_self.GQPoints[k][1] - 1.0                              #dN3/dl2
#   _self.dNdl2[k][3] = -4.0*_self.GQPoints[k][0]                                    #dN4/dl2
#   _self.dNdl2[k][4] =  4.0*_self.GQPoints[k][0]                                    #dN5/dl2
#   _self.dNdl2[k][5] =  4.0 - 4.0*_self.GQPoints[k][0] - 8.0*_self.GQPoints[k][1]   #dN6/dl2

   # Coordinate Transfomation
   # Lewis pag. 64 Eq. 3.108 for 1D
   # Lewis pag. 66 Eq. 3.121 for 2D
   for i in range(0,_self.NUMNODE):
    ii = _self.IEN[_e][i]
    
    dzdl1[k] += _self.z[ii]*_self.dNdl1[k][i]   # dx/dl1 = x1*dN1/dl1 + x2*dN2/dl1 + x3*dN3/dl1 ...
    drdl1[k] += _self.r[ii]*_self.dNdl1[k][i]   # dy/dl1 = y1*dN1/dl1 + y2*dN2/dl1 + y3*dN3/dl1 ...
    dzdl2[k] += _self.z[ii]*_self.dNdl2[k][i]   # dx/dl2 = x1*dN1/dl2 + x2*dN2/dl2 + x3*dN3/dl2 ...
    drdl2[k] += _self.r[ii]*_self.dNdl2[k][i]   # dy/dl2 = y1*dN1/dl2 + y2*dN2/dl2 + y3*dN3/dl2 ...


   # Jacobian Matrix
   # Lewis pag. 66 Eq. 3.121
   J[0][0] = dzdl1[k]
   J[0][1] = drdl1[k]
   J[1][0] = dzdl2[k]
   J[1][1] = drdl2[k]

   jacobian[k] = np.linalg.det(J)


   # Shape Functions Derivatives in respect to x and y
   # Lewis pag. 65 Eq. 3.116
   for i in range(0,_self.NUMNODE):
    dNdz[k][i] = (1.0/jacobian[k])*( _self.dNdl1[k][i]*drdl2[k] - _self.dNdl2[k][i]*drdl1[k])
    dNdr[k][i] = (1.0/jacobian[k])*(-_self.dNdl1[k][i]*dzdl2[k] + _self.dNdl2[k][i]*dzdl1[k])



  _self.mass = np.zeros([_self.NUMNODE,_self.NUMNODE], dtype = float)
  _self.kzz = np.zeros([_self.NUMNODE,_self.NUMNODE], dtype = float)
  _self.kzr = np.zeros([_self.NUMNODE,_self.NUMNODE], dtype = float)
  _self.krz = np.zeros([_self.NUMNODE,_self.NUMNODE], dtype = float)
  _self.krr = np.zeros([_self.NUMNODE,_self.NUMNODE], dtype = float)
  _self.gz = np.zeros([_self.NUMNODE,_self.NUMNODE], dtype = float)
  _self.gr = np.zeros([_self.NUMNODE,_self.NUMNODE], dtype = float)
  _self.dz = np.zeros([_self.NUMNODE,_self.NUMNODE], dtype = float)
  _self.dr = np.zeros([_self.NUMNODE,_self.NUMNODE], dtype = float)
  

  # Elementary Matrices
  # P.S: It is divided by 2 due to relation 
  # of parallelogram and triangle areas --> DxDy = (jacobian*Dl1*Dl2)/2
  for k in range(0,_self.NUMGAUSS): 
   for i in range(0,_self.NUMNODE):
    for j in range(0,_self.NUMNODE):
    
     _self.mass[i][j] += _self.N[k][i]*_self.N[k][j]*jacobian[k]*_self.GQWeights[k]/2.0
    
     _self.kzz[i][j] += dNdz[k][i]*dNdz[k][j]*jacobian[k]*_self.GQWeights[k]/2.0
     _self.kzr[i][j] += dNdz[k][i]*dNdr[k][j]*jacobian[k]*_self.GQWeights[k]/2.0
     _self.krz[i][j] += dNdr[k][i]*dNdz[k][j]*jacobian[k]*_self.GQWeights[k]/2.0
     _self.krr[i][j] += dNdr[k][i]*dNdr[k][j]*jacobian[k]*_self.GQWeights[k]/2.0
    
     _self.gz[i][j] += dNdz[k][j]*_self.N[k][i]*jacobian[k]*_self.GQWeights[k]/2.0
     _self.gr[i][j] += dNdr[k][j]*_self.N[k][i]*jacobian[k]*_self.GQWeights[k]/2.0
    
     _self.dz[i][j] += dNdz[k][j]*_self.N[k][i]*jacobian[k]*_self.GQWeights[k]/2.0
     _self.dr[i][j] += dNdr[k][j]*_self.N[k][i]*jacobian[k]*_self.GQWeights[k]/2.0
  



 def cubic(_self,_e):
  _self.NUMNODE = 10  #Cubic Triangular Element - 10 Nodes
  
  N = np.zeros([_self.NUMGAUSS,_self.NUMNODE], dtype = float)

  dNdl1 = np.zeros([_self.NUMGAUSS,_self.NUMNODE], dtype = float)
  dNdl2 = np.zeros([_self.NUMGAUSS,_self.NUMNODE], dtype = float)
  
  dzdl1 = np.zeros([_self.NUMGAUSS,1], dtype = float)
  dzdl2 = np.zeros([_self.NUMGAUSS,1], dtype = float)
  drdl1 = np.zeros([_self.NUMGAUSS,1], dtype = float)
  drdl2 = np.zeros([_self.NUMGAUSS,1], dtype = float)

  dNdz = np.zeros([_self.NUMGAUSS,_self.NUMNODE], dtype = float)
  dNdr = np.zeros([_self.NUMGAUSS,_self.NUMNODE], dtype = float)

  J = np.zeros([2,2], dtype = float)
  jacobian = np.zeros([_self.NUMGAUSS,1], dtype = float)




  # A loop is required for each pair of coordinates (l1,l2)
  for k in range(0,_self.NUMGAUSS):
    
   # Area Coordinates
   L1 = _self.GQPoints[k][0]                                #L1 = l1
   L2 = _self.GQPoints[k][1]                                #L2 = l2
   L3 = 1.0 - _self.GQPoints[k][0] - _self.GQPoints[k][1]   #L3 = 1 - l1 - l2

   # Shape Functions
   # Lewis pag. 56 and 57 Eq. 3.81/3.82/3.83/3.84/3.85
   N[k][0] = 0.5*L1*(3.0*L1 - 1.0)*(3.0*L1 - 2.0)   #N1  = 1/2*L1*(3*L1-1)*(3*L1-2)
   N[k][1] = 0.5*L2*(3.0*L2 - 1.0)*(3.0*L2 - 2.0)   #N2  = 1/2*L2*(3*L2-1)*(3*L2-2)
   N[k][2] = 0.5*L3*(3.0*L3 - 1.0)*(3.0*L3 - 2.0)   #N3  = 1/2*L3*(3*L3-1)*(3*L3-2)
   N[k][3] = 4.5*L1*L2*(3.0*L1 - 1.0)               #N4  = 9/2*L1*L2(3*L1-1)
   N[k][4] = 4.5*L1*L2*(3.0*L2 - 1.0)               #N5  = 9/2*L1*L2(3*L2-1)
   N[k][5] = 4.5*L2*L3*(3.0*L2 - 1.0)               #N6  = 9/2*L2*L3(3*L2-1)
   N[k][6] = 4.5*L2*L3*(3.0*L3 - 1.0)               #N7  = 9/2*L2*L3(3*L3-1)
   N[k][7] = 4.5*L3*L1*(3.0*L3 - 1.0)               #N8  = 9/2*L3*L1(3*L3-1)
   N[k][8] = 4.5*L3*L1*(3.0*L1 - 1.0)               #N9  = 9/2*L3*L1(3*L1-1)
   N[k][9] = 27.0*L1*L2*L3                          #N10 = 27*L1*L2*L3


   # Shape Functions Derivatives in respect to l1
   #dN1/dl1
   dNdl1[k][0] =  0.5*( 27.0*_self.GQPoints[k][0]**2 - 18.0*_self.GQPoints[k][0] + 2.0)

   #dN2/dl1
   dNdl1[k][1] =  0.0

   #dN3/dl1
   dNdl1[k][2] =  0.5*(-11.0 + 36.0*_self.GQPoints[k][1] - 54.0*_self.GQPoints[k][0]*_self.GQPoints[k][1]\
                      - 27.0*_self.GQPoints[k][1]**2 + 36.0*_self.GQPoints[k][0]\
                      - 27.0*_self.GQPoints[k][0]**2)
   #dN4/dl1
   dNdl1[k][3] = 4.5*( 6.0*_self.GQPoints[k][0]*_self.GQPoints[k][1] - _self.GQPoints[k][1])

   #dN5/dl1
   dNdl1[k][4] = 4.5*( 3.0*_self.GQPoints[k][1]**2 - _self.GQPoints[k][1])

   #dN6/dl1
   dNdl1[k][5] = 4.5*(-3.0*_self.GQPoints[k][1]**2 + _self.GQPoints[k][1])

   #dN7/dl1
   dNdl1[k][6] = 4.5*(-5.0*_self.GQPoints[k][1] + 6.0*_self.GQPoints[k][0]*_self.GQPoints[k][1]\
                     + 6.0*_self.GQPoints[k][1]**2)

   #dN8/dl1
   dNdl1[k][7] = 4.5*( 2.0 - 5.0*_self.GQPoints[k][1] - 10.0*_self.GQPoints[k][0]\
                     + 3.0*_self.GQPoints[k][1]**2 + 12.0*_self.GQPoints[k][0]*_self.GQPoints[k][1]\
                     + 9.0*_self.GQPoints[k][0]**2)

   #dN9/dl1
   dNdl1[k][8] = 4.5*(-1.0 + _self.GQPoints[k][1] + 8.0*_self.GQPoints[k][0]\
                     - 6.0*_self.GQPoints[k][0]*_self.GQPoints[k][1] - 9.0*_self.GQPoints[k][0]**2)

   #dN10/dl1
   dNdl1[k][9] = 27.0*(_self.GQPoints[k][1] - 2.0*_self.GQPoints[k][0]*_self.GQPoints[k][1]\
                     - _self.GQPoints[k][1]**2)



   # Shape Functions Derivatives in respect to l2
   #dN1/dl2
   dNdl2[k][0] =  0.0

   #dN2/dl2
   dNdl2[k][1] =  0.5*( 27.0*_self.GQPoints[k][1]**2 - 18.0*_self.GQPoints[k][1] + 2.0)

   #dN3/dl2
   dNdl2[k][2] =  0.5*(-11.0 + 36.0*_self.GQPoints[k][0] - 54.0*_self.GQPoints[k][0]*_self.GQPoints[k][1]\
                      - 27.0*_self.GQPoints[k][0]**2 + 36.0*_self.GQPoints[k][1]\
                      - 27.0*_self.GQPoints[k][1]**2)

   #dN4/dl2
   dNdl2[k][3] = 4.5*( 3.0*_self.GQPoints[k][0]**2 - _self.GQPoints[k][0])

   #dN5/dl2
   dNdl2[k][4] = 4.5*( 6.0*_self.GQPoints[k][0]*_self.GQPoints[k][1] - _self.GQPoints[k][0])

   #dN6/dl2
   dNdl2[k][5] = 4.5*(-1.0 + _self.GQPoints[k][0] + 8.0*_self.GQPoints[k][1]\
                     - 6.0*_self.GQPoints[k][0]*_self.GQPoints[k][1] - 9.0*_self.GQPoints[k][1]**2)

   #dN7/dl2
   dNdl2[k][6] = 4.5*( 2.0 - 5.0*_self.GQPoints[k][0] - 10.0*_self.GQPoints[k][1]\
                     + 3.0*_self.GQPoints[k][0]**2 + 12.0*_self.GQPoints[k][0]*_self.GQPoints[k][1]\
                     + 9.0*_self.GQPoints[k][1]**2)

   #dN8/dl2
   dNdl2[k][7] = 4.5*(-5.0*_self.GQPoints[k][0] + 6.0*_self.GQPoints[k][0]*_self.GQPoints[k][1]\
                     + 6.0*_self.GQPoints[k][0]**2)

   #dN9/dl2
   dNdl2[k][8] = 4.5*(-3.0*_self.GQPoints[k][0]**2 + _self.GQPoints[k][0])

   #dN10/dl2
   dNdl2[k][9] = 27.0*(_self.GQPoints[k][0] - 2.0*_self.GQPoints[k][0]*_self.GQPoints[k][1]\
                     - _self.GQPoints[k][0]**2)



   # Coordinate Transfomation
   # Lewis pag. 64 Eq. 3.108 for 1D
   # Lewis pag. 66 Eq. 3.121 for 2D
   for i in range(0,_self.NUMNODE):
    ii = _self.IEN[_e][i]
    
    dzdl1[k] += _self.z[ii]*dNdl1[k][i]   # dx/dl1 = x1*dN1/dl1 + x2*dN2/dl1 + x3*dN3/dl1 ...
    drdl1[k] += _self.r[ii]*dNdl1[k][i]   # dy/dl1 = y1*dN1/dl1 + y2*dN2/dl1 + y3*dN3/dl1 ...
    dzdl2[k] += _self.z[ii]*dNdl2[k][i]   # dx/dl2 = x1*dN1/dl2 + x2*dN2/dl2 + x3*dN3/dl2 ...
    drdl2[k] += _self.r[ii]*dNdl2[k][i]   # dy/dl2 = y1*dN1/dl2 + y2*dN2/dl2 + y3*dN3/dl2 ...



   # Jacobian Matrix
   # Lewis pag. 66 Eq. 3.121
   J[0][0] = dzdl1[k]
   J[0][1] = drdl1[k]
   J[1][0] = dzdl2[k]
   J[1][1] = drdl2[k]

   jacobian[k] = np.linalg.det(J)


   # Shape Functions Derivatives in respect to x and y
   # Lewis pag. 65 Eq. 3.116
   for i in range(0,_self.NUMNODE):
    dNdz[k][i] = (1.0/jacobian[k])*( dNdl1[k][i]*drdl2[k] - dNdl2[k][i]*drdl1[k])
    dNdr[k][i] = (1.0/jacobian[k])*(-dNdl1[k][i]*dzdl2[k] + dNdl2[k][i]*dzdl1[k])



  _self.mass = np.zeros([_self.NUMNODE,_self.NUMNODE], dtype = float)
  _self.kzz = np.zeros([_self.NUMNODE,_self.NUMNODE], dtype = float)
  _self.kzr = np.zeros([_self.NUMNODE,_self.NUMNODE], dtype = float)
  _self.krz = np.zeros([_self.NUMNODE,_self.NUMNODE], dtype = float)
  _self.krr = np.zeros([_self.NUMNODE,_self.NUMNODE], dtype = float)
  _self.gz = np.zeros([_self.NUMNODE,_self.NUMNODE], dtype = float)
  _self.gr = np.zeros([_self.NUMNODE,_self.NUMNODE], dtype = float)
  _self.dz = np.zeros([_self.NUMNODE,_self.NUMNODE], dtype = float)
  _self.dr = np.zeros([_self.NUMNODE,_self.NUMNODE], dtype = float)
  

  # Elementary Matrices
  # P.S: It is divided by 2 due to relation 
  # of parallelogram and triangle areas --> DxDy = (jacobian*Dl1*Dl2)/2
  for k in range(0,_self.NUMGAUSS): 
   for i in range(0,_self.NUMNODE):
    for j in range(0,_self.NUMNODE):
    
     _self.mass[i][j] += N[k][i]*N[k][j]*jacobian[k]*_self.GQWeights[k]/2.0
    
     _self.kzz[i][j] += dNdz[k][i]*dNdz[k][j]*jacobian[k]*_self.GQWeights[k]/2.0
     _self.kzr[i][j] += dNdz[k][i]*dNdr[k][j]*jacobian[k]*_self.GQWeights[k]/2.0
     _self.krz[i][j] += dNdr[k][i]*dNdz[k][j]*jacobian[k]*_self.GQWeights[k]/2.0
     _self.krr[i][j] += dNdr[k][i]*dNdr[k][j]*jacobian[k]*_self.GQWeights[k]/2.0
    
     _self.gz[i][j] += dNdz[k][j]*N[k][i]*jacobian[k]*_self.GQWeights[k]/2.0
     _self.gr[i][j] += dNdr[k][j]*N[k][i]*jacobian[k]*_self.GQWeights[k]/2.0
    
     _self.dz[i][j] += dNdz[k][j]*N[k][i]*jacobian[k]*_self.GQWeights[k]/2.0
     _self.dr[i][j] += dNdr[k][j]*N[k][i]*jacobian[k]*_self.GQWeights[k]/2.0



 def assembleMv(_self,_e):
  _self.NUMNODE = 3  #Linear Triangular Element - 3 Nodes
  
  N = np.zeros([_self.NUMGAUSS,_self.NUMNODE], dtype = float)

  dNdl1 = np.zeros([_self.NUMGAUSS,_self.NUMNODE], dtype = float)
  dNdl2 = np.zeros([_self.NUMGAUSS,_self.NUMNODE], dtype = float)
  
  dzdl1 = np.zeros([_self.NUMGAUSS,1], dtype = float)
  dzdl2 = np.zeros([_self.NUMGAUSS,1], dtype = float)
  drdl1 = np.zeros([_self.NUMGAUSS,1], dtype = float)
  drdl2 = np.zeros([_self.NUMGAUSS,1], dtype = float)

  dNdz = np.zeros([_self.NUMGAUSS,_self.NUMNODE], dtype = float)
  dNdr = np.zeros([_self.NUMGAUSS,_self.NUMNODE], dtype = float)

  J = np.zeros([2,2], dtype = float)
  jacobian = np.zeros([_self.NUMGAUSS,1], dtype = float)


  # A loop is required for each pair of coordinates (l1,l2)
  for k in range(0,_self.NUMGAUSS):
    
   # Area Coordinates
   # Lewis pag. 67 Eq. 3.129
   L1 = 1.0 - _self.GQPoints[k][0] - _self.GQPoints[k][1]   #L1 = 1 - l1 - l2
   L2 = _self.GQPoints[k][0]                                #L2 = l1
   L3 = _self.GQPoints[k][1]                                #L3 = l2

   # Shape Functions
   # Lewis pag. 67 Eq. 3.129
   N[k][0] = L1  #N1 = L1
   N[k][1] = L2  #N2 = L2
   N[k][2] = L3  #N3 = L3

   # Shape Functions Derivatives in respect to l1
   dNdl1[k][0] = -1.0   #dN1/dl1
   dNdl1[k][1] =  1.0   #dN2/dl1
   dNdl1[k][2] =  0.0   #dN3/dl1

   # Shape Functions Derivatives in respect to l2
   dNdl2[k][0] = -1.0   #dN1/dl2
   dNdl2[k][1] =  0.0   #dN2/dl2
   dNdl2[k][2] =  1.0   #dN3/dl2

   # Coordinate Transfomation
   # Lewis pag. 64 Eq. 3.108 for 1D
   # Lewis pag. 66 Eq. 3.121 for 2D
   for i in range(0,_self.NUMNODE):
    ii = _self.IEN[_e][i]
    
    dzdl1[k] += _self.z[ii]*dNdl1[k][i]   # dz/dl1 = z1*dN1/dl1 + z2*dN2/dl1 + z3*dN3/dl1 ...
    drdl1[k] += _self.r[ii]*dNdl1[k][i]   # dr/dl1 = r1*dN1/dl1 + r2*dN2/dl1 + r3*dN3/dl1 ...
    dzdl2[k] += _self.z[ii]*dNdl2[k][i]   # dz/dl2 = z1*dN1/dl2 + z2*dN2/dl2 + z3*dN3/dl2 ...
    drdl2[k] += _self.r[ii]*dNdl2[k][i]   # dr/dl2 = r1*dN1/dl2 + r2*dN2/dl2 + r3*dN3/dl2 ...


   # Jacobian Matrix
   # Lewis pag. 66 Eq. 3.121
   J[0][0] = dzdl1[k]
   J[0][1] = drdl1[k]
   J[1][0] = dzdl2[k]
   J[1][1] = drdl2[k]

   jacobian[k] = np.linalg.det(J)


   # Shape Functions Derivatives in respect to x and y
   # Lewis pag. 65 Eq. 3.116
   for i in range(0,_self.NUMNODE):
    dNdz[k][i] = (1.0/jacobian[k])*( dNdl1[k][i]*drdl2[k] - dNdl2[k][i]*drdl1[k])
    dNdr[k][i] = (1.0/jacobian[k])*(-dNdl1[k][i]*dzdl2[k] + dNdl2[k][i]*dzdl1[k])



  _self.mass = np.zeros([_self.NUMNODE,_self.NUMNODE], dtype = float)
  
  # Elementary Matrices
  # P.S: It is divided by 2 due to relation 
  # of parallelogram and triangle areas --> DxDy = (jacobian*Dl1*Dl2)/2
  for k in range(0,_self.NUMGAUSS): 
   for i in range(0,_self.NUMNODE):
    for j in range(0,_self.NUMNODE):
    
     _self.mass[i][j] += N[k][i]*N[k][j]*jacobian[k]*_self.GQWeights[k]/2.0







 # Reference Analytic: Fundamentals of the Finite
 #                     Element Method for Heat Transfer
 #                     and Fluid Flow - Lewis, Nithiarasu,
 #                     Seetharamu - pg. 196-200
 #                     For Q_elem pg. 126
 #                     For 1D pg. 193
 def analytic2D(_self, _e):

  i = _self.IEN[_e][0]
  j = _self.IEN[_e][1]
  k = _self.IEN[_e][2]

  bi = _self.r[j]-_self.r[k]
  bj = _self.r[k]-_self.r[i]
  bk = _self.r[i]-_self.r[j]
  ci = _self.z[k]-_self.z[j]
  cj = _self.z[i]-_self.z[k]
  ck = _self.z[j]-_self.z[i]


  A = 0.5*np.linalg.det(np.array([[1, _self.z[i], _self.r[i]],
 				  [1, _self.z[j], _self.r[j]],
				  [1, _self.z[k], _self.r[k]]]))


  _self.mass = (A/12.)*np.array([[2.,1.,1.],
                                 [1.,2.,1.],
                                 [1.,1.,2.]])

  _self.q = (A/3.)*np.ones([3,1], dtype = float)
  
  _self.gz = (1./6)*np.array([[bi,bj,bk],
                              [bi,bj,bk],
                              [bi,bj,bk]]) 
   
  _self.gr = (1./6)*np.array([[ci,cj,ck],
                              [ci,cj,ck],
                              [ci,cj,ck]])

  _self.kzz = (1./(4*A))*np.array([[bi*bi,bi*bj,bi*bk],
                                   [bj*bi,bj*bj,bj*bk],
                                   [bk*bi,bk*bj,bk*bk]])

  _self.krr = (1./(4*A))*np.array([[ci*ci,ci*cj,ci*ck],
                                   [cj*ci,cj*cj,cj*ck],
                                   [ck*ci,ck*cj,ck*ck]])

  _self.kzr = (1./(4*A))*np.array([[bi*ci,bi*cj,bi*ck],
                                   [bj*ci,bj*cj,bj*ck],
                                   [bk*ci,bk*cj,bk*ck]])

  _self.krz = (1./(4*A))*np.array([[ci*bi,ci*bj,ci*bk],
                                   [cj*bi,cj*bj,cj*bk],
                                   [ck*bi,ck*bj,ck*bk]])



 def analytic_axi(_self, _e):
  _self.r = _self.r
  _self.z = _self.z

  i = _self.IEN[_e][0]
  j = _self.IEN[_e][1]
  k = _self.IEN[_e][2]

  bi = _self.z[j] - _self.z[k]
  bj = _self.z[k] - _self.z[i]
  bk = _self.z[i] - _self.z[j]
  ci = _self.r[k] - _self.r[j]
  cj = _self.r[i] - _self.r[k]
  ck = _self.r[j] - _self.r[i]

  A = 0.5*np.linalg.det(np.array([[1, _self.r[i], _self.z[i]],
 				  [1, _self.r[j], _self.z[j]],
				  [1, _self.r[k], _self.z[k]]]))

  r = (_self.r[i] + _self.r[j] + _self.r[k])/3.

  r_vec = np.array([[_self.r[i]],
                    [_self.r[j]],
                    [_self.r[k]]])

  _self.M_elem = (A/12.)*np.array([[2.,1.,1.],
				   [1.,2.,1.],
				   [1.,1.,2.]])

  _self.Q_elem = (2*np.pi)*np.dot(_self.M_elem,r_vec)
  
  _self.Gr_elem = (1./6)*np.array([[bi,bj,bk],
                                   [bi,bj,bk],
                                   [bi,bj,bk]]) 
   
  _self.Gz_elem = (1./6)*np.array([[ci,cj,ck],
                                   [ci,cj,ck],
                                   [ci,cj,ck]])

  _self.Kr_elem = ((2*np.pi*r)/(4*A))*np.array([[bi*bi,bj*bi,bk*bi],
                                                [bi*bj,bj*bj,bk*bj],
                                                [bi*bk,bj*bk,bk*bk]])

  _self.Kz_elem = ((2*np.pi*r)/(4*A))*np.array([[ci*ci,cj*ci,ck*ci],
                                                [ci*cj,cj*cj,ck*cj],
                                                [ci*ck,cj*ck,ck*ck]])






