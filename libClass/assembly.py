# ==========================================
# Code created by Leandro Marques at 12/2018
# Gesar Search Group
# State University of the Rio de Janeiro
# e-mail: marquesleandro67@gmail.com
# ==========================================

# This code is used to assembly global matrices

# ------------------------------------------------------------------------------
# Use:
# Kxx, Kxy, Kyx, Kyy, K, M, MLump, Gx, Gy = 
# assembly.Linear2D(mesh.GL, mesh.npoints, mesh.nelem, mesh.IEN, mesh.x, mesh.y)
# ------------------------------------------------------------------------------

import sys
import numpy as np
import gaussianQuadrature
import scipy.sparse as sps
from tqdm import tqdm



def Element1D(_polynomial_option, _GL, _npoints, _nelem, _IEN, _x, _GAUSSPOINTS):
 K = sps.lil_matrix((_npoints,_npoints), dtype = float)
 M = sps.lil_matrix((_npoints,_npoints), dtype = float)
 G = sps.lil_matrix((_npoints,_npoints), dtype = float)
 
 element1D = gaussianQuadrature.Element1D(_x, _IEN, _GAUSSPOINTS)
 
 if _polynomial_option == 1:
  polynomial_order = 'Linear Element'
  
  for e in tqdm(range(0, _nelem)):
   element1D.linear(e)

   for i in range(0,_GL): 
    ii = _IEN[e][i]
  
    for j in range(0,_GL):
     jj = _IEN[e][j]

     K[ii,jj] += element1D.kx[i][j]
     M[ii,jj] += element1D.mass[i][j]
     G[ii,jj] += element1D.gx[i][j]


 elif _polynomial_option == 2:
  polynomial_order = 'Quadratic Element'

  for e in tqdm(range(0, _nelem)):
   element1D.quadratic(e)

   for i in range(0,_GL): 
    ii = _IEN[e][i]
  
    for j in range(0,_GL):
     jj = _IEN[e][j]

     K[ii,jj] += element1D.kx[i][j]
     M[ii,jj] += element1D.mass[i][j]
     G[ii,jj] += element1D.gx[i][j]



 else:
  print ""
  print " Error: Element type not found"
  print ""
  sys.exit()


 return K, M, G, polynomial_order

 



def Element2D(_polynomial_option, _GL, _npoints, _nelem, _IEN, _x, _y, _GAUSSPOINTS):

 Kxx = sps.lil_matrix((_npoints,_npoints), dtype = float)
 Kxy = sps.lil_matrix((_npoints,_npoints), dtype = float)
 Kyx = sps.lil_matrix((_npoints,_npoints), dtype = float)
 Kyy = sps.lil_matrix((_npoints,_npoints), dtype = float)
 K = sps.lil_matrix((_npoints,_npoints), dtype = float)
 M = sps.lil_matrix((_npoints,_npoints), dtype = float)
 MLump = sps.lil_matrix((_npoints,_npoints), dtype = float)
 Gx = sps.lil_matrix((_npoints,_npoints), dtype = float)
 Gy = sps.lil_matrix((_npoints,_npoints), dtype = float)


 element2D = gaussianQuadrature.Element2D(_x, _y, _IEN, _GAUSSPOINTS)

 if _polynomial_option == 1:
  polynomial_order = 'Linear Element'
  
  for e in tqdm(range(0, _nelem)):
   element2D.linear(e)

   for i in range(0,_GL): 
    ii = _IEN[e][i]
  
    for j in range(0,_GL):
     jj = _IEN[e][j]

     Kxx[ii,jj] += element2D.kxx[i][j]
     Kxy[ii,jj] += element2D.kxy[i][j]
     Kyx[ii,jj] += element2D.kyx[i][j]
     Kyy[ii,jj] += element2D.kyy[i][j]
     K[ii,jj] += element2D.kxx[i][j] + element2D.kyy[i][j]
   
     M[ii,jj] += element2D.mass[i][j]
     MLump[ii,ii] += element2D.mass[i][j]

     Gx[ii,jj] += element2D.gx[i][j]
     Gy[ii,jj] += element2D.gy[i][j]

 elif _polynomial_option == 2:
  polynomial_order = 'Mini Element'

  for e in tqdm(range(0, _nelem)):
   element2D.mini(e)

   for i in range(0,_GL): 
    ii = _IEN[e][i]
  
    for j in range(0,_GL):
     jj = _IEN[e][j]

     Kxx[ii,jj] += element2D.kxx[i][j]
     Kxy[ii,jj] += element2D.kxy[i][j]
     Kyx[ii,jj] += element2D.kyx[i][j]
     Kyy[ii,jj] += element2D.kyy[i][j]
     K[ii,jj] += element2D.kxx[i][j] + element2D.kyy[i][j]
   
     M[ii,jj] += element2D.mass[i][j]
     MLump[ii,ii] += element2D.mass[i][j]

     Gx[ii,jj] += element2D.gx[i][j]
     Gy[ii,jj] += element2D.gy[i][j]


 elif _polynomial_option == 3:
  polynomial_order = 'Quadratic Element'

  for e in tqdm(range(0, _nelem)):
   element2D.quadratic(e)

   for i in range(0,_GL): 
    ii = _IEN[e][i]
  
    for j in range(0,_GL):
     jj = _IEN[e][j]

     Kxx[ii,jj] += element2D.kxx[i][j]
     Kxy[ii,jj] += element2D.kxy[i][j]
     Kyx[ii,jj] += element2D.kyx[i][j]
     Kyy[ii,jj] += element2D.kyy[i][j]
     K[ii,jj] += element2D.kxx[i][j] + element2D.kyy[i][j]
   
     M[ii,jj] += element2D.mass[i][j]
     MLump[ii,ii] += element2D.mass[i][j]

     Gx[ii,jj] += element2D.gx[i][j]
     Gy[ii,jj] += element2D.gy[i][j]

 elif _polynomial_option == 4:
  polynomial_order = 'Cubic Element'

  for e in tqdm(range(0, _nelem)):
   element2D.cubic(e)

   for i in range(0,_GL): 
    ii = _IEN[e][i]
  
    for j in range(0,_GL):
     jj = _IEN[e][j]

     Kxx[ii,jj] += element2D.kxx[i][j]
     Kxy[ii,jj] += element2D.kxy[i][j]
     Kyx[ii,jj] += element2D.kyx[i][j]
     Kyy[ii,jj] += element2D.kyy[i][j]
     K[ii,jj] += element2D.kxx[i][j] + element2D.kyy[i][j]
   
     M[ii,jj] += element2D.mass[i][j]
     MLump[ii,ii] += element2D.mass[i][j]

     Gx[ii,jj] += element2D.gx[i][j]
     Gy[ii,jj] += element2D.gy[i][j]



 else:
  print ""
  print " Error: Element type not found"
  print ""
  sys.exit()


 return Kxx, Kxy, Kyx, Kyy, K, M, MLump, Gx, Gy, polynomial_order


def AxiElement2D(_simulation_option, _polynomial_option, _GL, _npoints, _nelem, _IEN, _z, _r, _GAUSSPOINTS):

 Kzzr = sps.lil_matrix((_npoints,_npoints), dtype = float)
 Krrr = sps.lil_matrix((_npoints,_npoints), dtype = float)
 Mr = sps.lil_matrix((_npoints,_npoints), dtype = float)
 M1r = sps.lil_matrix((_npoints,_npoints), dtype = float)
 Mr2 = sps.lil_matrix((_npoints,_npoints), dtype = float)
 Gr   = sps.lil_matrix((_npoints,_npoints), dtype = float)
 Gz   = sps.lil_matrix((_npoints,_npoints), dtype = float)
 Grr   = sps.lil_matrix((_npoints,_npoints), dtype = float)
 Gzr   = sps.lil_matrix((_npoints,_npoints), dtype = float)

 element2D = gaussianQuadrature.AxiElement2D(_z, _r, _IEN, _GAUSSPOINTS)
 
 if _simulation_option == 1:
  if _polynomial_option == 1:
   polynomial_order = 'Linear Element'
   
   for e in tqdm(range(0, _nelem)):
    element2D.linear(e)
 
    v1 = _IEN[e][0]
    v2 = _IEN[e][1]
    v3 = _IEN[e][2]
 
    r_elem = (_r[v1] + _r[v2] + _r[v3])/3.0
 
 
    for i in range(0,_GL): 
     ii = _IEN[e][i]
   
     for j in range(0,_GL):
      jj = _IEN[e][j]
 
      Kzzr[ii,jj] += r_elem*element2D.kzz[i][j]
      Krrr[ii,jj] += r_elem*element2D.krr[i][j]
    
      Mr[ii,jj] += r_elem*element2D.mass[i][j]
      M1r[ii,jj] += (1.0/r_elem)*element2D.mass[i][j]
      Mr2[ii,jj] += (r_elem**2)*element2D.mass[i][j]
 
      Gr[ii,jj] += element2D.gr[i][j]
      Gz[ii,jj] += element2D.gz[i][j]
      Grr[ii,jj] += r_elem*element2D.gr[i][j]
      Gzr[ii,jj] += r_elem*element2D.gz[i][j]
 
 
 
  elif _polynomial_option == 2:
   polynomial_order = 'Mini Element'
   
   for e in tqdm(range(0, _nelem)):
    element2D.mini(e)
 
    v1 = _IEN[e][0]
    v2 = _IEN[e][1]
    v3 = _IEN[e][2]
    v4 = _IEN[e][3]
 
    r_elem = (_r[v1] + _r[v2] + _r[v3] + _r[v4])/4.0
 
 
    for i in range(0,_GL): 
     ii = _IEN[e][i]
   
     for j in range(0,_GL):
      jj = _IEN[e][j]
 
      Kzzr[ii,jj] += r_elem*element2D.kzz[i][j]
      Krrr[ii,jj] += r_elem*element2D.krr[i][j]
    
      Mr[ii,jj] += r_elem*element2D.mass[i][j]
      M1r[ii,jj] += (1.0/r_elem)*element2D.mass[i][j]
      Mr2[ii,jj] += (r_elem**2)*element2D.mass[i][j]
 
      Gr[ii,jj] += element2D.gr[i][j]
      Gz[ii,jj] += element2D.gz[i][j]
      Grr[ii,jj] += r_elem*element2D.gr[i][j]
      Gzr[ii,jj] += r_elem*element2D.gz[i][j]
 
 
 
  elif _polynomial_option == 3:
   polynomial_order = 'Quad Element'
   
   for e in tqdm(range(0, _nelem)):
    element2D.quadratic(e)
 
    v1 = _IEN[e][0]
    v2 = _IEN[e][1]
    v3 = _IEN[e][2]
    v4 = _IEN[e][3]
    v5 = _IEN[e][4]
    v6 = _IEN[e][5]
 
    r_elem = (_r[v1] + _r[v2] + _r[v3] + _r[v4] + _r[v5] + _r[v6])/6.0
   
    for i in range(0,_GL): 
     ii = _IEN[e][i]
   
     for j in range(0,_GL):
      jj = _IEN[e][j]
 
      Kzzr[ii,jj] += r_elem*element2D.kzz[i][j]
      Krrr[ii,jj] += r_elem*element2D.krr[i][j]
    
      Mr[ii,jj] += r_elem*element2D.mass[i][j]
      M1r[ii,jj] += (1.0/r_elem)*element2D.mass[i][j]
      Mr2[ii,jj] += (r_elem**2)*element2D.mass[i][j]
 
      Gr[ii,jj] += element2D.gr[i][j]
      Gz[ii,jj] += element2D.gz[i][j]
      Grr[ii,jj] += r_elem*element2D.gr[i][j]
      Gzr[ii,jj] += r_elem*element2D.gz[i][j]
 
 
  else:
   print ""
   print " Error: Element type not found"
   print ""
   sys.exit()

 
 #Debug
 elif _simulation_option == 2:
  polynomial_order = 'Debug'
  Kzzr = Kzzr*1.0  
  Krrr = Krrr*1.0
  Mr   = Mr*1.0
  M1r  = M1r*1.0
  Mr2  = Mr2*1.0
  Gr   = Gr*1.0
  Gz   = Gz*1.0
  Grr  = Grr*1.0  
  Gzr  = Gzr*1.0
 


 return Kzzr, Krrr, Mr, M1r, Mr2, Gr, Gz, Grr, Gzr, polynomial_order


def AxiAssembleMv(_polynomial_option, _GL, _npoints, _nelem, _IEN, _z, _r, _vr, _GAUSSPOINTS):

 Mv = sps.lil_matrix((_npoints,_npoints), dtype = float)

 element2D = gaussianQuadrature.AxiElement2D(_z, _r, _IEN, _GAUSSPOINTS)
 
 if _polynomial_option == 1:
  polynomial_order = 'Linear Element'
  
  for e in tqdm(range(0, _nelem)):
   element2D.assembleMv_linear(e)

   v1 = _IEN[e][0]
   v2 = _IEN[e][1]
   v3 = _IEN[e][2]

   v_ele = (_vr[v1] + _vr[v2] + _vr[v3])/3.0


   for i in range(0,_GL): 
    ii = _IEN[e][i]
  
    for j in range(0,_GL):
     jj = _IEN[e][j]

     Mv[ii,jj] += v_ele*element2D.mass[i][j]



 elif _polynomial_option == 2:
  polynomial_order = 'Mini Element'
  
  for e in tqdm(range(0, _nelem)):
   element2D.assembleMv_mini(e)

   v1 = _IEN[e][0]
   v2 = _IEN[e][1]
   v3 = _IEN[e][2]
   v4 = _IEN[e][3]

   v_ele = (_vr[v1] + _vr[v2] + _vr[v3] + _vr[v4])/4.0


   for i in range(0,_GL): 
    ii = _IEN[e][i]
  
    for j in range(0,_GL):
     jj = _IEN[e][j]

     Mv[ii,jj] += v_ele*element2D.mass[i][j]


 elif _polynomial_option == 3:
  polynomial_order = 'Quad Element'
  
  for e in tqdm(range(0, _nelem)):
   element2D.assembleMv_quadratic(e)

   v1 = _IEN[e][0]
   v2 = _IEN[e][1]
   v3 = _IEN[e][2]
   v4 = _IEN[e][3]
   v5 = _IEN[e][4]
   v6 = _IEN[e][5]

   v_ele = (_vr[v1] + _vr[v2] + _vr[v3] + _vr[v4] + _vr[v5] + _vr[v6])/6.0
   #v_ele = (_vr[v1] + _vr[v2] + _vr[v3])/3.0
  
   for i in range(0,_GL): 
    ii = _IEN[e][i]
  
    for j in range(0,_GL):
     jj = _IEN[e][j]

     Mv[ii,jj] += v_ele*element2D.mass[i][j]

 else:
  print ""
  print " Error: Element type not found"
  print ""
  sys.exit()


 return Mv






def AxiElement2D_domega(_polynomial_option, _GL, _npoints, _nelem, _IEN, _z, _r, _GAUSSPOINTS):

 Kzz = sps.lil_matrix((_npoints,_npoints), dtype = float)
 Kzr = sps.lil_matrix((_npoints,_npoints), dtype = float)
 Krz = sps.lil_matrix((_npoints,_npoints), dtype = float)
 Krr = sps.lil_matrix((_npoints,_npoints), dtype = float)
 K = sps.lil_matrix((_npoints,_npoints), dtype = float)
 M = sps.lil_matrix((_npoints,_npoints), dtype = float)
 Mr = sps.lil_matrix((_npoints,_npoints), dtype = float)
 M1r = sps.lil_matrix((_npoints,_npoints), dtype = float)
 M1r2 = sps.lil_matrix((_npoints,_npoints), dtype = float)
 MLump = sps.lil_matrix((_npoints,_npoints), dtype = float)
 Gz   = sps.lil_matrix((_npoints,_npoints), dtype = float)
 Gz1r = sps.lil_matrix((_npoints,_npoints), dtype = float)
 Gr   = sps.lil_matrix((_npoints,_npoints), dtype = float)
 Gr1r = sps.lil_matrix((_npoints,_npoints), dtype = float)


 element2D = gaussianQuadrature.Element2D(_z, _r, _IEN, _GAUSSPOINTS)
 
 if _polynomial_option == 1:
  polynomial_order = 'Linear Element'
  
  for e in tqdm(range(0, _nelem)):
   element2D.linear(e)

   v1 = _IEN[e][0]
   v2 = _IEN[e][1]
   v3 = _IEN[e][2]

   r_ele = (_r[v1] + _r[v2] + _r[v3])/3.0


   for i in range(0,_GL): 
    ii = _IEN[e][i]
  
    for j in range(0,_GL):
     jj = _IEN[e][j]

     Kzz[ii,jj] += element2D.kxx[i][j]
     Kzr[ii,jj] += element2D.kxy[i][j]
     Krz[ii,jj] += element2D.kyx[i][j]
     Krr[ii,jj] += element2D.kyy[i][j]
     K[ii,jj] += element2D.kxx[i][j] + element2D.kyy[i][j]
   
     M[ii,jj] += element2D.mass[i][j]
     Mr[ii,jj] += r_ele*element2D.mass[i][j]
     M1r[ii,jj] += (1.0/r_ele)*element2D.mass[i][j]
     M1r2[ii,jj] += (1.0/(r_ele**2))*element2D.mass[i][j]
     MLump[ii,ii] += element2D.mass[i][j]

     Gz[ii,jj] += element2D.gx[i][j]
     Gr[ii,jj] += element2D.gy[i][j]
     Gz1r[ii,jj] += (1.0/r_ele)*element2D.gx[i][j]
     Gr1r[ii,jj] += (1.0/r_ele)*element2D.gy[i][j]



 elif _polynomial_option == 2:
  polynomial_order = 'Mini Element'
  
  for e in tqdm(range(0, _nelem)):
   element2D.mini(e)

   v1 = _IEN[e][0]
   v2 = _IEN[e][1]
   v3 = _IEN[e][2]
   v4 = _IEN[e][3]

   r_ele = (_r[v1] + _r[v2] + _r[v3])/3.0


   for i in range(0,_GL): 
    ii = _IEN[e][i]
  
    for j in range(0,_GL):
     jj = _IEN[e][j]

     Kzz[ii,jj] += element2D.kxx[i][j]
     Kzr[ii,jj] += element2D.kxy[i][j]
     Krz[ii,jj] += element2D.kyx[i][j]
     Krr[ii,jj] += element2D.kyy[i][j]
     K[ii,jj] += element2D.kxx[i][j] + element2D.kyy[i][j]
   
     M[ii,jj] += element2D.mass[i][j]
     Mr[ii,jj] += r_ele*element2D.mass[i][j]
     M1r[ii,jj] += (1.0/r_ele)*element2D.mass[i][j]
     M1r2[ii,jj] += (1.0/(r_ele**2))*element2D.mass[i][j]
     MLump[ii,ii] += element2D.mass[i][j]

     Gz[ii,jj] += element2D.gx[i][j]
     Gr[ii,jj] += element2D.gy[i][j]
     Gz1r[ii,jj] += (1.0/r_ele)*element2D.gx[i][j]
     Gr1r[ii,jj] += (1.0/r_ele)*element2D.gy[i][j]





 elif _polynomial_option == 3:
  polynomial_order = 'Quad Element'
  
  for e in tqdm(range(0, _nelem)):
   element2D.quadratic(e)

   v1 = _IEN[e][0]
   v2 = _IEN[e][1]
   v3 = _IEN[e][2]
   v4 = _IEN[e][3]
   v5 = _IEN[e][4]
   v6 = _IEN[e][5]

   #r_ele = (_r[v1] + _r[v2] + _r[v3] + _r[v4] + _r[v5] + _r[v6])/6.0
   r_ele = (_r[v1] + _r[v2] + _r[v3])/3.0
  
   for i in range(0,_GL): 
    ii = _IEN[e][i]
  
    for j in range(0,_GL):
     jj = _IEN[e][j]

     Kzz[ii,jj] += element2D.kxx[i][j]
     Kzr[ii,jj] += element2D.kxy[i][j]
     Krz[ii,jj] += element2D.kyx[i][j]
     Krr[ii,jj] += element2D.kyy[i][j]
     K[ii,jj] += element2D.kxx[i][j] + element2D.kyy[i][j]
   
     M[ii,jj] += element2D.mass[i][j]
     Mr[ii,jj] += r_ele*element2D.mass[i][j]
     M1r[ii,jj] += (1.0/r_ele)*element2D.mass[i][j]
     M1r2[ii,jj] += (1.0/(r_ele**2))*element2D.mass[i][j]
     MLump[ii,ii] += element2D.mass[i][j]

     Gz[ii,jj] += element2D.gx[i][j]
     Gr[ii,jj] += element2D.gy[i][j]
     Gz1r[ii,jj] += (1.0/r_ele)*element2D.gx[i][j]
     Gr1r[ii,jj] += (1.0/r_ele)*element2D.gy[i][j]

 else:
  print ""
  print " Error: Element type not found"
  print ""
  sys.exit()


 return Kzz, Kzr, Krz, Krr, K, M, Mr, M1r, M1r2, MLump, Gz, Gr, Gz1r, Gr1r, polynomial_order




def Mini_NS2D(_GLV, _GLP, _NV, _NP, _nelem, _IEN, _x, _y):
 
 K = sps.lil_matrix((2*_NV,2*_NV), dtype = float)
 M = sps.lil_matrix((2*_NV,2*_NV), dtype = float)
 MLump = sps.lil_matrix((2*_NV,2*_NV), dtype = float)
 G = sps.lil_matrix((2*_NV,_NP), dtype = float)
 D = sps.lil_matrix((_NP,2*_NV), dtype = float)


 mini = gaussianQuadrature.Mini(_x, _y, _IEN)

 for e in tqdm(range(0, _nelem)):
  mini.numerical(e)

  for i in range(0, _GLV): 
   ii = _IEN[e][i]
  
   for j in range(0, _GLV):
    jj = _IEN[e][j]

    #MSC 2007 pag.84
    K[ii,jj] += 2.0*mini.kxx[i][j] + mini.kyy[i][j] #K11
    K[ii,jj + _NV] += mini.kxy[i][j] #K12
    K[ii + _NV,jj] += mini.kyx[i][j] #K21
    K[ii + _NV,jj + _NV] += mini.kxx[i][j] + 2.0*mini.kyy[i][j] #K22
   
    M[ii,jj] += mini.mass[i][j]
    M[ii + _NV,jj + _NV] += mini.mass[i][j]
    
    MLump[ii,ii] += mini.mass[i][j]
    MLump[ii + _NV,ii + _NV] += mini.mass[i][j]


   for k in range(0, _GLP):
    kk = _IEN[e][k]

    G[ii,kk] += mini.gx[i][k]
    G[ii + _NV,kk] += mini.gy[i][k]

    D[kk,ii] += mini.dx[k][i]
    D[kk,ii + _NV] += mini.dy[k][i]


 return K, M, MLump, G, D



