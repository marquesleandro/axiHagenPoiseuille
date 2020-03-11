# ==========================================
# Code created by Leandro Marques at 03/2020
# Gesar Search Group
# State University of the Rio de Janeiro
# e-mail: marquesleandro67@gmail.com
# ==========================================

# This code is used to import .vtk file


# Converting .msh in a python list


def vtkfile(_file,_npoints): 

 vtklist = [] 
 with open(_file) as vtkfile:
   for line in vtkfile:
    row = line.split()
    vtklist.append(row[:])

 if vtklist.find("scalar1"):
  print vtklist.index("scalar1")
  print "yes"
 
 '''
 for i in range(0,_npoints):
  if vtklist[i][0] == "POINTS":
   line_points = i + 1

  elif vtklist[i][0] == "VECTORS":
   line_vector = i + 1

  elif vtklist[i][1] == "scalar1":
   line_scalar1 = i + 1

  elif vtklist[i][1] == "scalar2":
   line_scalar2 = i + 1

  elif vtklist[i][1] == "scalar3":
   line_scalar3 = i + 1


 x = np.zeros([_npoints,1], dtype = float)
 y = np.zeros([_npoints,1], dtype = float)
 vx = np.zeros([_npoints,1], dtype = float)
 vy = np.zeros([_npoints,1], dtype = float)
 scalar1 = np.zeros([_npoints,1], dtype = float)
 scalar2 = np.zeros([_npoints,1], dtype = float)
 scalar3 = np.zeros([_npoints,1], dtype = float)
 for i in range(0,len(_npoints)):
  x[i] = vtk[line_points + i][0]
  y[i] = vtk[line_points + i][1]

  
  vx[i] = vtk[line_vector + i][0]
  vy[i] = vtk[line_vector + i][1]

  scalar1[i] = vtk[line_scalar1 + i][0]
  scalar2[i] = vtk[line_scalar2 + i][0]
  scalar3[i] = vtk[line_scalar3 + i][0]

 return x, y, vx, vy, scalar1, scalar2, scalar3
 '''
