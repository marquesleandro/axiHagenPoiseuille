# ==========================================
# Code created by Leandro Marques at 12/2018
# Gesar Search Group
# State University of the Rio de Janeiro
# e-mail: marquesleandro67@gmail.com
# ==========================================

# This code is used to import .msh
# for n governament equations 


import numpy as np


class Linear1D:
 def __init__(_self, _dir, _file, _number_equations):
  _self.name = _dir + '/' + _file
  _self.neq = _number_equations
  _self.mshFile = []
  _self.neumann_lines = {}
  _self.neumann_pts = {}
  _self.dirichlet_lines = {}
  _self.dirichlet_pts = {}
  _self.neighborsNodes = {}
  _self.neighborsElements = {}
  _self.far_neighborsNodes = {}
  _self.far_neighborsElements = {}

  # Converting .msh in a python list
  with open(_self.name) as mesh:
   for line in mesh:
    row = line.split()
    _self.mshFile.append(row[:])

  
  _self.numPhysical = int(_self.mshFile[4][0])
  _self.numNodes = int(_self.mshFile[_self.numPhysical+7][0])
  
  for i in range(1, _self.neq + 1):
   _self.neumann_lines[i] = []
   _self.neumann_pts[i] = []
   _self.dirichlet_lines[i] = []
   _self.dirichlet_pts[i] = []

  # Lines classification in neumann or dirichlet
  for i in range(0,_self.numPhysical):
   for j in range(1,_self.neq + 1):
    if _self.mshFile[5+i][2] == '"neumann' + str(j):
     _self.neumann_lines[j].append(int(_self.mshFile[5+i][1]))

    elif _self.mshFile[5+i][2] == '"dirichlet' + str(j):
     _self.dirichlet_lines[j].append(int(_self.mshFile[5+i][1]))


  countLineStart = _self.numPhysical + _self.numNodes 
  countLine = 1
  mshFileElementType = int(_self.mshFile[(countLineStart + 10 + countLine)][1])

  # Assembly of the neumann edges and dirichlet points
  while mshFileElementType == 15:
    a_1 = int(_self.mshFile[(countLineStart + 10 + countLine)][3])
    a_2 = int(_self.mshFile[(countLineStart + 10 + countLine)][5])
    a_3 = [a_1,a_2]
   
    for i in range(1, _self.neq + 1): 
     if a_1 in _self.neumann_lines[i]:
      _self.neumann_pts[i].append(a_3)
 
      countLine = countLine + 1
      mshFileElementType = int(_self.mshFile[(countLineStart + 10 + countLine)][1])
 
     elif a_1 in _self.dirichlet_lines[i]:
      _self.dirichlet_pts[i].append(a_3)

      countLine = countLine + 1
      mshFileElementType = int(_self.mshFile[(countLineStart + 10 + countLine)][1])

   
  for i in range(1, _self.neq + 1):
   _self.neumann_pts[i] = np.array(_self.neumann_pts[i]) 
   _self.dirichlet_pts[i] = np.array(_self.dirichlet_pts[i]) 


  _self.numElements = int(_self.mshFile[countLineStart + 10][0]) - countLine + 1

  for i in range(0, _self.numNodes):
   _self.neighborsNodes[i] = []
   _self.neighborsElements[i] = []
   _self.far_neighborsNodes[i] = []
   _self.far_neighborsElements[i] = []


  countLine = countLine
  countLineStart = countLineStart


 def coord(_self):
  _self.x = np.zeros([_self.numNodes,1], dtype = float)
  _self.npts = []

  for i in range(0, _self.numNodes):  
   _self.x[i] = _self.mshFile[_self.numPhysical + 8 + i][1]
   _self.npts.append(i)


 def ien(_self):
  _self.IEN = np.zeros([_self.numElements,2], dtype = int)
  _self.FreedomDegree = len(_self.IEN[0,:])
  length = [] 

  for e in range(0, _self.numElements):
   v1 = int(_self.mshFile[(countLineStart + 10 + countLine + e)][5]) - 1
   v2 = int(_self.mshFile[(countLineStart + 10 + countLine + e)][6]) - 1
  
   _self.IEN[e] = [v1,v2]
  
   _self.neighborsNodes[v1].extend(_self.IEN[e])  
   _self.neighborsNodes[v2].extend(_self.IEN[e])  
   
   _self.neighborsNodes[v1] = list(set(_self.neighborsNodes[v1]))
   _self.neighborsNodes[v2] = list(set(_self.neighborsNodes[v2]))
   
   _self.neighborsElements[v1].append(e)  
   _self.neighborsElements[v2].append(e)  

   x_a = _self.x[v1] - _self.x[v2]
   length1 = np.sqrt(x_a**2)
   length.append(length1)

  _self.minLengthMesh = min(length)


  for i in range(0, _self.numNodes):
   for j in _self.neighborsNodes[i]:
    _self.far_neighborsNodes[i].extend(_self.neighborsNodes[j]) 
    _self.far_neighborsElements[i].extend(_self.neighborsElements[j]) 
 
   _self.far_neighborsNodes[i] = list(set(_self.far_neighborsNodes[i])\
                                     - set(_self.neighborsNodes[i]))
   
   _self.far_neighborsElements[i] = list(set(_self.far_neighborsElements[i])\
                                        - set(_self.neighborsElements[i]))



class Quad1D:
 def __init__(_self, _dir, _file, _number_equations):
  _self.name = _dir + '/' + _file
  _self.neq = _number_equations
  _self.mshFile = []
  _self.neumann_lines = {}
  _self.neumann_pts = {}
  _self.dirichlet_lines = {}
  _self.dirichlet_pts = {}
  _self.neighborsNodes = {}
  _self.neighborsElements = {}
  _self.far_neighborsNodes = {}
  _self.far_neighborsElements = {}

  # Converting .msh in a python list
  with open(_self.name) as mesh:
   for line in mesh:
    row = line.split()
    _self.mshFile.append(row[:])

  
  _self.numPhysical = int(_self.mshFile[4][0])
  _self.numNodes = int(_self.mshFile[_self.numPhysical+7][0])
  
  for i in range(1, _self.neq + 1):
   _self.neumann_lines[i] = []
   _self.neumann_pts[i] = []
   _self.dirichlet_lines[i] = []
   _self.dirichlet_pts[i] = []

  # Lines classification in neumann or dirichlet
  for i in range(0,_self.numPhysical):
   for j in range(1,_self.neq + 1):
    if _self.mshFile[5+i][2] == '"neumann' + str(j):
     _self.neumann_lines[j].append(int(_self.mshFile[5+i][1]))

    elif _self.mshFile[5+i][2] == '"dirichlet' + str(j):
     _self.dirichlet_lines[j].append(int(_self.mshFile[5+i][1]))


  countLineStart = _self.numPhysical + _self.numNodes 
  countLine = 1
  mshFileElementType = int(_self.mshFile[(countLineStart + 10 + countLine)][1])

  # Assembly of the neumann edges and dirichlet points
  while mshFileElementType == 15:
    a_1 = int(_self.mshFile[(countLineStart + 10 + countLine)][3])
    a_2 = int(_self.mshFile[(countLineStart + 10 + countLine)][5])
    a_3 = [a_1,a_2]
   
    for i in range(1, _self.neq + 1): 
     if a_1 in _self.neumann_lines[i]:
      _self.neumann_pts[i].append(a_3)
 
      countLine = countLine + 1
      mshFileElementType = int(_self.mshFile[(countLineStart + 10 + countLine)][1])
 
     elif a_1 in _self.dirichlet_lines[i]:
      _self.dirichlet_pts[i].append(a_3)

      countLine = countLine + 1
      mshFileElementType = int(_self.mshFile[(countLineStart + 10 + countLine)][1])

   
  for i in range(1, _self.neq + 1):
   _self.neumann_pts[i] = np.array(_self.neumann_pts[i]) 
   _self.dirichlet_pts[i] = np.array(_self.dirichlet_pts[i]) 


  _self.numElements = int(_self.mshFile[countLineStart + 10][0]) - countLine + 1

  for i in range(0, _self.numNodes):
   _self.neighborsNodes[i] = []
   _self.neighborsElements[i] = []
   _self.far_neighborsNodes[i] = []
   _self.far_neighborsElements[i] = []


  countLine = countLine
  countLineStart = countLineStart


 def coord(_self):
  _self.x = np.zeros([_self.numNodes,1], dtype = float)
  _self.npts = []

  for i in range(0, _self.numNodes):  
   _self.x[i] = _self.mshFile[_self.numPhysical + 8 + i][1]
   _self.npts.append(i)


 def ien(_self):
  _self.IEN = np.zeros([_self.numElements,3], dtype = int)
  _self.FreedomDegree = len(_self.IEN[0,:])
  length = [] 

  for e in range(0, _self.numElements):
   v1 = int(_self.mshFile[(countLineStart + 10 + countLine + e)][5]) - 1
   v2 = int(_self.mshFile[(countLineStart + 10 + countLine + e)][6]) - 1
   v3 = int(_self.mshFile[(countLineStart + 10 + countLine + e)][7]) - 1
  
   _self.IEN[e] = [v1,v2,v3]
  
   _self.neighborsNodes[v1].extend(_self.IEN[e])  
   _self.neighborsNodes[v2].extend(_self.IEN[e])  
   _self.neighborsNodes[v3].extend(_self.IEN[e])  
   
   _self.neighborsNodes[v1] = list(set(_self.neighborsNodes[v1]))
   _self.neighborsNodes[v2] = list(set(_self.neighborsNodes[v2]))
   _self.neighborsNodes[v3] = list(set(_self.neighborsNodes[v2]))

   _self.neighborsElements[v1].append(e)  
   _self.neighborsElements[v2].append(e)  
   _self.neighborsElements[v3].append(e)  

   x_a = _self.x[v1] - _self.x[v2]
   length1 = np.sqrt(x_a**2)
   length.append(length1)

  _self.minLengthMesh = min(length)


  for i in range(0, _self.numNodes):
   for j in _self.neighborsNodes[i]:
    _self.far_neighborsNodes[i].extend(_self.neighborsNodes[j]) 
    _self.far_neighborsElements[i].extend(_self.neighborsElements[j]) 
 
   _self.far_neighborsNodes[i] = list(set(_self.far_neighborsNodes[i])\
                                     - set(_self.neighborsNodes[i]))
   
   _self.far_neighborsElements[i] = list(set(_self.far_neighborsElements[i])\
                                        - set(_self.neighborsElements[i]))




class Linear2D:
 def __init__(_self, _dir, _file):
  _self.name = _dir + '/' + _file
  _self.mshFile = []
  _self.boundaryEdges = []
  _self.neighborsNodes = {}
  _self.neighborsElements = {}


  
  # Converting .msh in a python list
  with open(_self.name) as mesh:
   for line in mesh:
    row = line.split()
    _self.mshFile.append(row[:])

  _self.numPhysical = int(_self.mshFile[4][0])
  _self.numNodes = int(_self.mshFile[_self.numPhysical+7][0])
 
 

  # Boundary edges Assembly
  countLineStart = _self.numPhysical + _self.numNodes 
  countLine = 1
  mshFileElementType = int(_self.mshFile[(countLineStart + 10 + countLine)][1])

  while mshFileElementType == 1:
   a_1 = int(_self.mshFile[(countLineStart + 10 + countLine)][3])
   a_2 = int(_self.mshFile[(countLineStart + 10 + countLine)][5])
   a_3 = int(_self.mshFile[(countLineStart + 10 + countLine)][6])
   a_4 = [a_1,a_2,a_3]
   
   _self.boundaryEdges.append(a_4)
 
   countLine = countLine + 1
   mshFileElementType = int(_self.mshFile[(countLineStart + 10 + countLine)][1])
 
  _self.boundaryEdges = np.array(_self.boundaryEdges) 
  _self.numElements = int(_self.mshFile[countLineStart + 10][0]) - countLine + 1



  # Coordinate vectors Assembly
  _self.x = np.zeros([_self.numNodes,1], dtype = float)
  _self.y = np.zeros([_self.numNodes,1], dtype = float)
  _self.npts = []

  for i in range(0, _self.numNodes):  
   _self.x[i] = _self.mshFile[_self.numPhysical + 8 + i][1]
   _self.y[i] = _self.mshFile[_self.numPhysical + 8 + i][2]
   _self.neighborsNodes[i] = []
   _self.neighborsElements[i] = []
   _self.npts.append(i)



  # IEN matrix and Neighbors list Assembly
  _self.IEN = np.zeros([_self.numElements,3], dtype = int)
  _self.FreedomDegree = len(_self.IEN[0,:])
  length = []

  for e in range(0, _self.numElements):
   v1 = int(_self.mshFile[(countLineStart + 10 + countLine + e)][5]) - 1
   v2 = int(_self.mshFile[(countLineStart + 10 + countLine + e)][6]) - 1
   v3 = int(_self.mshFile[(countLineStart + 10 + countLine + e)][7]) - 1
  
   _self.IEN[e] = [v1,v2,v3]
  
   _self.neighborsNodes[v1].extend(_self.IEN[e])  
   _self.neighborsNodes[v2].extend(_self.IEN[e])  
   _self.neighborsNodes[v3].extend(_self.IEN[e])  
   
   _self.neighborsNodes[v1] = list(set(_self.neighborsNodes[v1]))
   _self.neighborsNodes[v2] = list(set(_self.neighborsNodes[v2]))
   _self.neighborsNodes[v3] = list(set(_self.neighborsNodes[v3]))
   
   _self.neighborsElements[v1].append(e)  
   _self.neighborsElements[v2].append(e)  
   _self.neighborsElements[v3].append(e)  

   x_a = _self.x[v1] - _self.x[v2]
   x_b = _self.x[v2] - _self.x[v3]
   x_c = _self.x[v3] - _self.x[v1]
   
   y_a = _self.y[v1] - _self.y[v2]
   y_b = _self.y[v2] - _self.y[v3]
   y_c = _self.y[v3] - _self.y[v1]
   
   length1 = np.sqrt(x_a**2 + y_a**2)
   length2 = np.sqrt(x_b**2 + y_b**2)
   length3 = np.sqrt(x_c**2 + y_c**2)

   length.append(length1)
   length.append(length2)
   length.append(length3)
   
  _self.minLengthMesh = min(length)



class Mini2D:
 def __init__(_self, _dir, _file, _number_equations):
  _self.name = _dir + '/' + _file
  _self.neq = _number_equations
  _self.mshFile = []
  _self.neumann_lines = {}
  _self.dirichlet_lines = {}
  _self.neumann_edges = {}
  _self.dirichlet_pts = {}
  _self.neighborsNodes = {}
  _self.neighborsElements = {}
  _self.neighborsNodes_linear = {}
  _self.neighborsElements_linear = {}
  _self.far_neighborsNodes = {}
  _self.far_neighborsElements = {}

  
  # Converting .msh in a python list
  with open(_self.name) as mesh:
   for line in mesh:
    row = line.split()
    _self.mshFile.append(row[:])

  
  _self.numPhysical = int(_self.mshFile[4][0])
  _self.numNodes = int(_self.mshFile[_self.numPhysical+7][0])
  
  for i in range(1, _self.neq + 1):
   _self.neumann_lines[i] = []
   _self.neumann_edges[i] = []
   _self.dirichlet_lines[i] = []
   _self.dirichlet_pts[i] = []


  # Lines classification in neumann or dirichlet
  for i in range(0,_self.numPhysical):
   for j in range(1,_self.neq + 1):
    if _self.mshFile[5+i][2] == '"neumann' + str(j):
     _self.neumann_lines[j].append(int(_self.mshFile[5+i][1]))

    elif _self.mshFile[5+i][2] == '"dirichlet' + str(j):
     _self.dirichlet_lines[j].append(int(_self.mshFile[5+i][1]))


  countLineStart = _self.numPhysical + _self.numNodes 
  countLine = 1
  mshFileElementType = int(_self.mshFile[(countLineStart + 10 + countLine)][1])

  # Assembly of the neumann edges and dirichlet points
  while mshFileElementType == 1:
    a_1 = int(_self.mshFile[(countLineStart + 10 + countLine)][3])
    a_2 = int(_self.mshFile[(countLineStart + 10 + countLine)][5])
    a_3 = int(_self.mshFile[(countLineStart + 10 + countLine)][6])
    a_4 = [a_1,a_2,a_3]
   
    for i in range(1, _self.neq + 1): 
     if a_1 in _self.neumann_lines[i]:
      _self.neumann_edges[i].append(a_4)
 
      countLine = countLine + 1
      mshFileElementType = int(_self.mshFile[(countLineStart + 10 + countLine)][1])
 
     elif a_1 in _self.dirichlet_lines[i]:
      _self.dirichlet_pts[i].append(a_4)

      countLine = countLine + 1
      mshFileElementType = int(_self.mshFile[(countLineStart + 10 + countLine)][1])

   
  for i in range(1, _self.neq + 1):
   _self.neumann_edges[i] = np.array(_self.neumann_edges[i]) 
   _self.dirichlet_pts[i] = np.array(_self.dirichlet_pts[i]) 


  _self.numElements = int(_self.mshFile[countLineStart + 10][0]) - countLine + 1
  _self.numNodes_linear = _self.numNodes
  _self.numNodes = _self.numNodes_linear + _self.numElements


  for i in range(0, _self.numNodes):
   _self.neighborsNodes[i] = []
   _self.neighborsElements[i] = []
   _self.far_neighborsNodes[i] = []
   _self.far_neighborsElements[i] = []

  for i in range(0, _self.numNodes_linear):
   _self.neighborsNodes_linear[i] = []
   _self.neighborsElements_linear[i] = []
 

  countLine = countLine
  countLineStart = countLineStart


 def coord(_self):
  _self.x = np.zeros([_self.numNodes,1], dtype = float)
  _self.y = np.zeros([_self.numNodes,1], dtype = float)
  _self.x_linear = np.zeros([_self.numNodes_linear,1], dtype = float)
  _self.y_linear = np.zeros([_self.numNodes_linear,1], dtype = float)
  _self.npts = []

  for i in range(0, _self.numNodes_linear):  
   _self.x[i] = _self.mshFile[_self.numPhysical + 8 + i][1]
   _self.y[i] = _self.mshFile[_self.numPhysical + 8 + i][2]
   _self.x_linear[i] = _self.mshFile[_self.numPhysical + 8 + i][1]
   _self.y_linear[i] = _self.mshFile[_self.numPhysical + 8 + i][2]
   _self.npts.append(i)


 def ien(_self):
  _self.IEN = np.zeros([_self.numElements,4], dtype = int)
  _self.IEN_linear = np.zeros([_self.numElements,3], dtype = int)
  _self.FreedomDegree = len(_self.IEN[0,:])
  _self.FreedomDegree_linear = len(_self.IEN_linear[0,:])
  length = []
  _self.nodes_linear = []

  for e in range(0, _self.numElements):
   v1 = int(_self.mshFile[(countLineStart + 10 + countLine + e)][5]) - 1
   v2 = int(_self.mshFile[(countLineStart + 10 + countLine + e)][6]) - 1
   v3 = int(_self.mshFile[(countLineStart + 10 + countLine + e)][7]) - 1
   v4 = _self.numNodes_linear + e
  
   _self.IEN[e] = [v1,v2,v3,v4]
   _self.IEN_linear[e] = [v1,v2,v3]

   _self.x[v4] = (_self.x[v1] + _self.x[v2] + _self.x[v3])/3.0
   _self.y[v4] = (_self.y[v1] + _self.y[v2] + _self.y[v3])/3.0
   _self.npts.append(v4)

   _self.nodes_linear.append(v1)  
   _self.nodes_linear.append(v2)  
   _self.nodes_linear.append(v3)  

   _self.neighborsNodes[v1].extend(_self.IEN[e])  
   _self.neighborsNodes[v2].extend(_self.IEN[e])  
   _self.neighborsNodes[v3].extend(_self.IEN[e])  
   _self.neighborsNodes[v4].extend(_self.IEN[e])  
   
   _self.neighborsNodes[v1] = list(set(_self.neighborsNodes[v1]))
   _self.neighborsNodes[v2] = list(set(_self.neighborsNodes[v2]))
   _self.neighborsNodes[v3] = list(set(_self.neighborsNodes[v3]))
   _self.neighborsNodes[v4] = list(set(_self.neighborsNodes[v4]))
   
   _self.neighborsElements[v1].append(e)  
   _self.neighborsElements[v2].append(e)  
   _self.neighborsElements[v3].append(e)  
   _self.neighborsElements[v4].append(e)  

   _self.neighborsNodes_linear[v1].extend(_self.IEN_linear[e])  
   _self.neighborsNodes_linear[v2].extend(_self.IEN_linear[e])  
   _self.neighborsNodes_linear[v3].extend(_self.IEN_linear[e])  
   
   _self.neighborsNodes_linear[v1] = list(set(_self.neighborsNodes_linear[v1]))
   _self.neighborsNodes_linear[v2] = list(set(_self.neighborsNodes_linear[v2]))
   _self.neighborsNodes_linear[v3] = list(set(_self.neighborsNodes_linear[v3]))
   
   _self.neighborsElements_linear[v1].append(e)  
   _self.neighborsElements_linear[v2].append(e)  
   _self.neighborsElements_linear[v3].append(e)  

   x_a = _self.x[v1] - _self.x[v2]
   x_b = _self.x[v2] - _self.x[v3]
   x_c = _self.x[v3] - _self.x[v1]
   
   y_a = _self.y[v1] - _self.y[v2]
   y_b = _self.y[v2] - _self.y[v3]
   y_c = _self.y[v3] - _self.y[v1]
   
   length1 = np.sqrt(x_a**2 + y_a**2)
   length2 = np.sqrt(x_b**2 + y_b**2)
   length3 = np.sqrt(x_c**2 + y_c**2)

   length.append(length1)
   length.append(length2)
   length.append(length3)
   
  _self.minLengthMesh = min(length)
  _self.nodes_linear = list(set(_self.nodes_linear))  

  for i in range(0, _self.numNodes):
   for j in _self.neighborsNodes[i]:
    _self.far_neighborsNodes[i].extend(_self.neighborsNodes[j]) 
    _self.far_neighborsElements[i].extend(_self.neighborsElements[j]) 
 
   _self.far_neighborsNodes[i] = list(set(_self.far_neighborsNodes[i])\
                                     - set(_self.neighborsNodes[i]))
   
   _self.far_neighborsElements[i] = list(set(_self.far_neighborsElements[i])\
                                        - set(_self.neighborsElements[i]))





class Quad2D:
 def __init__(_self, _dir, _file, _number_equations):
  _self.name = _dir + '/' + _file
  _self.neq = _number_equations
  _self.mshFile = []
  _self.neumann_lines = {}
  _self.neumann_edges = {}
  _self.dirichlet_lines = {}
  _self.dirichlet_pts = {}
  _self.neighborsNodes = {}
  _self.neighborsElements = {}
  _self.far_neighborsNodes = {}
  _self.far_neighborsElements = {}

  # Converting .msh in a python list
  with open(_self.name) as mesh:
   for line in mesh:
    row = line.split()
    _self.mshFile.append(row[:])

  
  _self.numPhysical = int(_self.mshFile[4][0])
  _self.numNodes = int(_self.mshFile[_self.numPhysical+7][0])
  
  for i in range(1, _self.neq + 1):
   _self.neumann_lines[i] = []
   _self.neumann_edges[i] = []
   _self.dirichlet_lines[i] = []
   _self.dirichlet_pts[i] = []


  # Lines classification in neumann or dirichlet
  for i in range(0,_self.numPhysical):
   for j in range(1,_self.neq + 1):
    if _self.mshFile[5+i][2] == '"neumann' + str(j):
     _self.neumann_lines[j].append(int(_self.mshFile[5+i][1]))

    elif _self.mshFile[5+i][2] == '"dirichlet' + str(j):
     _self.dirichlet_lines[j].append(int(_self.mshFile[5+i][1]))


  countLineStart = _self.numPhysical + _self.numNodes 
  countLine = 1
  mshFileElementType = int(_self.mshFile[(countLineStart + 10 + countLine)][1])

  # Assembly of the neumann edges and dirichlet points
  while mshFileElementType == 8:
    a_1 = int(_self.mshFile[(countLineStart + 10 + countLine)][3])
    a_2 = int(_self.mshFile[(countLineStart + 10 + countLine)][5])
    a_3 = int(_self.mshFile[(countLineStart + 10 + countLine)][6])
    a_4 = int(_self.mshFile[(countLineStart + 10 + countLine)][7])
    a_5 = [a_1,a_2,a_3,a_4]
   
    for i in range(1, _self.neq + 1): 
     if a_1 in _self.neumann_lines[i]:
      _self.neumann_edges[i].append(a_5)
 
      countLine = countLine + 1
      mshFileElementType = int(_self.mshFile[(countLineStart + 10 + countLine)][1])
 
     elif a_1 in _self.dirichlet_lines[i]:
      _self.dirichlet_pts[i].append(a_5)

      countLine = countLine + 1
      mshFileElementType = int(_self.mshFile[(countLineStart + 10 + countLine)][1])

   
  for i in range(1, _self.neq + 1):
   _self.neumann_edges[i] = np.array(_self.neumann_edges[i]) 
   _self.dirichlet_pts[i] = np.array(_self.dirichlet_pts[i]) 


  _self.numElements = int(_self.mshFile[countLineStart + 10][0]) - countLine + 1

  for i in range(0, _self.numNodes):
   _self.neighborsNodes[i] = []
   _self.neighborsElements[i] = []
   _self.far_neighborsNodes[i] = []
   _self.far_neighborsElements[i] = []


  countLine = countLine
  countLineStart = countLineStart


 def coord(_self):
  _self.x = np.zeros([_self.numNodes,1], dtype = float)
  _self.y = np.zeros([_self.numNodes,1], dtype = float)
  _self.npts = []

  for i in range(0, _self.numNodes):  
   _self.x[i] = _self.mshFile[_self.numPhysical + 8 + i][1]
   _self.y[i] = _self.mshFile[_self.numPhysical + 8 + i][2]
   _self.npts.append(i)


 def ien(_self):
  _self.IEN = np.zeros([_self.numElements,6], dtype = int)
  _self.FreedomDegree = len(_self.IEN[0,:])
  _self.nodes_linear = [] 
  _self.nodes_quad = [] 
  length = [] 
  
  for e in range(0, _self.numElements):
   v1 = int(_self.mshFile[(countLineStart + 10 + countLine + e)][5]) - 1
   v2 = int(_self.mshFile[(countLineStart + 10 + countLine + e)][6]) - 1
   v3 = int(_self.mshFile[(countLineStart + 10 + countLine + e)][7]) - 1
   v4 = int(_self.mshFile[(countLineStart + 10 + countLine + e)][8]) - 1
   v5 = int(_self.mshFile[(countLineStart + 10 + countLine + e)][9]) - 1
   v6 = int(_self.mshFile[(countLineStart + 10 + countLine + e)][10]) - 1
  
   _self.IEN[e] = [v1,v2,v3,v4,v5,v6]
  
   _self.neighborsNodes[v1].extend(_self.IEN[e])  
   _self.neighborsNodes[v2].extend(_self.IEN[e])  
   _self.neighborsNodes[v3].extend(_self.IEN[e])  
   _self.neighborsNodes[v4].extend(_self.IEN[e])  
   _self.neighborsNodes[v5].extend(_self.IEN[e])  
   _self.neighborsNodes[v6].extend(_self.IEN[e])  
   
   _self.neighborsNodes[v1] = list(set(_self.neighborsNodes[v1]))
   _self.neighborsNodes[v2] = list(set(_self.neighborsNodes[v2]))
   _self.neighborsNodes[v3] = list(set(_self.neighborsNodes[v3]))
   _self.neighborsNodes[v4] = list(set(_self.neighborsNodes[v4]))
   _self.neighborsNodes[v5] = list(set(_self.neighborsNodes[v5]))
   _self.neighborsNodes[v6] = list(set(_self.neighborsNodes[v6]))
   
   _self.neighborsElements[v1].append(e)  
   _self.neighborsElements[v2].append(e)  
   _self.neighborsElements[v3].append(e)  
   _self.neighborsElements[v4].append(e)  
   _self.neighborsElements[v5].append(e)  
   _self.neighborsElements[v6].append(e)  
   
   _self.nodes_linear.append(v1)  
   _self.nodes_linear.append(v2)  
   _self.nodes_linear.append(v3)  

   _self.nodes_quad.append(v4)  
   _self.nodes_quad.append(v5)  
   _self.nodes_quad.append(v6)  

   x_a = _self.x[v1] - _self.x[v4]
   x_b = _self.x[v4] - _self.x[v2]
   x_c = _self.x[v2] - _self.x[v5]
   x_d = _self.x[v5] - _self.x[v3]
   x_e = _self.x[v3] - _self.x[v6]
   x_f = _self.x[v6] - _self.x[v1]
   
   y_a = _self.y[v1] - _self.y[v4]
   y_b = _self.y[v4] - _self.y[v2]
   y_c = _self.y[v2] - _self.y[v5]
   y_d = _self.y[v5] - _self.y[v3]
   y_e = _self.y[v3] - _self.y[v6]
   y_f = _self.y[v6] - _self.y[v1]
   
   length1 = np.sqrt(x_a**2 + y_a**2)
   length2 = np.sqrt(x_b**2 + y_b**2)
   length3 = np.sqrt(x_c**2 + y_c**2)
   length4 = np.sqrt(x_d**2 + y_d**2)
   length5 = np.sqrt(x_e**2 + y_e**2)
   length6 = np.sqrt(x_f**2 + y_f**2)

   length.append(length1)
   length.append(length2)
   length.append(length3)
   length.append(length4)
   length.append(length5)
   length.append(length6)
   
  _self.minLengthMesh = min(length)

  _self.nodes_linear = list(set(_self.nodes_linear))  
  _self.nodes_quad = list(set(_self.nodes_quad))


  for i in range(0, _self.numNodes):
   for j in _self.neighborsNodes[i]:
    _self.far_neighborsNodes[i].extend(_self.neighborsNodes[j]) 
    _self.far_neighborsElements[i].extend(_self.neighborsElements[j]) 
 
   _self.far_neighborsNodes[i] = list(set(_self.far_neighborsNodes[i])\
                                     - set(_self.neighborsNodes[i]))
   
   _self.far_neighborsElements[i] = list(set(_self.far_neighborsElements[i])\
                                        - set(_self.neighborsElements[i]))

   
class Cubic2D:
 def __init__(_self, _dir, _file, _number_equations):
  _self.name = _dir + '/' + _file
  _self.neq = _number_equations
  _self.mshFile = []
  _self.neumann_lines = {}
  _self.neumann_edges = {}
  _self.dirichlet_lines = {}
  _self.dirichlet_pts = {}
  _self.neighborsNodes = {}
  _self.neighborsElements = {}
  _self.far_neighborsNodes = {}
  _self.far_neighborsElements = {}

  # Converting .msh in a python list
  with open(_self.name) as mesh:
   for line in mesh:
    row = line.split()
    _self.mshFile.append(row[:])

  
  _self.numPhysical = int(_self.mshFile[4][0])
  _self.numNodes = int(_self.mshFile[_self.numPhysical+7][0])
  
  for i in range(1, _self.neq + 1):
   _self.neumann_lines[i] = []
   _self.neumann_edges[i] = []
   _self.dirichlet_lines[i] = []
   _self.dirichlet_pts[i] = []


  # Lines classification in neumann or dirichlet
  for i in range(0,_self.numPhysical):
   for j in range(1,_self.neq + 1):
    if _self.mshFile[5+i][2] == '"neumann' + str(j):
     _self.neumann_lines[j].append(int(_self.mshFile[5+i][1]))

    elif _self.mshFile[5+i][2] == '"dirichlet' + str(j):
     _self.dirichlet_lines[j].append(int(_self.mshFile[5+i][1]))


  countLineStart = _self.numPhysical + _self.numNodes 
  countLine = 1
  mshFileElementType = int(_self.mshFile[(countLineStart + 10 + countLine)][1])

  # Assembly of the neumann edges and dirichlet points
  while mshFileElementType == 26:
    a_1 = int(_self.mshFile[(countLineStart + 10 + countLine)][3])
    a_2 = int(_self.mshFile[(countLineStart + 10 + countLine)][5])
    a_3 = int(_self.mshFile[(countLineStart + 10 + countLine)][6])
    a_4 = int(_self.mshFile[(countLineStart + 10 + countLine)][7])
    a_5 = int(_self.mshFile[(countLineStart + 10 + countLine)][8])
    a_6 = [a_1,a_2,a_3,a_4,a_5]
   
    for i in range(1, _self.neq + 1): 
     if a_1 in _self.neumann_lines[i]:
      _self.neumann_edges[i].append(a_6)
 
      countLine = countLine + 1
      mshFileElementType = int(_self.mshFile[(countLineStart + 10 + countLine)][1])
 
     elif a_1 in _self.dirichlet_lines[i]:
      _self.dirichlet_pts[i].append(a_6)

      countLine = countLine + 1
      mshFileElementType = int(_self.mshFile[(countLineStart + 10 + countLine)][1])

   
  for i in range(1, _self.neq + 1):
   _self.neumann_edges[i] = np.array(_self.neumann_edges[i]) 
   _self.dirichlet_pts[i] = np.array(_self.dirichlet_pts[i]) 


  _self.numElements = int(_self.mshFile[countLineStart + 10][0]) - countLine + 1

  for i in range(0, _self.numNodes):
   _self.neighborsNodes[i] = []
   _self.neighborsElements[i] = []
   _self.far_neighborsNodes[i] = []
   _self.far_neighborsElements[i] = []


  countLine = countLine
  countLineStart = countLineStart


 def coord(_self):
  _self.x = np.zeros([_self.numNodes,1], dtype = float)
  _self.y = np.zeros([_self.numNodes,1], dtype = float)
  _self.npts = []

  for i in range(0, _self.numNodes):  
   _self.x[i] = _self.mshFile[_self.numPhysical + 8 + i][1]
   _self.y[i] = _self.mshFile[_self.numPhysical + 8 + i][2]
   _self.npts.append(i)


 def ien(_self):
  _self.IEN = np.zeros([_self.numElements,10], dtype = int)
  _self.FreedomDegree = len(_self.IEN[0,:])
  _self.nodes_linear = [] 
  _self.nodes_quad = [] 
  length = [] 
  
  for e in range(0, _self.numElements):
   v1 = int(_self.mshFile[(countLineStart + 10 + countLine + e)][5]) - 1
   v2 = int(_self.mshFile[(countLineStart + 10 + countLine + e)][6]) - 1
   v3 = int(_self.mshFile[(countLineStart + 10 + countLine + e)][7]) - 1
   v4 = int(_self.mshFile[(countLineStart + 10 + countLine + e)][8]) - 1
   v5 = int(_self.mshFile[(countLineStart + 10 + countLine + e)][9]) - 1
   v6 = int(_self.mshFile[(countLineStart + 10 + countLine + e)][10]) - 1
   v7 = int(_self.mshFile[(countLineStart + 10 + countLine + e)][11]) - 1
   v8 = int(_self.mshFile[(countLineStart + 10 + countLine + e)][12]) - 1
   v9 = int(_self.mshFile[(countLineStart + 10 + countLine + e)][13]) - 1
   v10 = int(_self.mshFile[(countLineStart + 10 + countLine + e)][14]) - 1
  
   _self.IEN[e] = [v1,v2,v3,v4,v5,v6,v7,v8,v9,v10]
  
   _self.neighborsNodes[v1].extend(_self.IEN[e])  
   _self.neighborsNodes[v2].extend(_self.IEN[e])  
   _self.neighborsNodes[v3].extend(_self.IEN[e])  
   _self.neighborsNodes[v4].extend(_self.IEN[e])  
   _self.neighborsNodes[v5].extend(_self.IEN[e])  
   _self.neighborsNodes[v6].extend(_self.IEN[e])  
   _self.neighborsNodes[v7].extend(_self.IEN[e])  
   _self.neighborsNodes[v8].extend(_self.IEN[e])  
   _self.neighborsNodes[v9].extend(_self.IEN[e])  
   _self.neighborsNodes[v10].extend(_self.IEN[e])  
   
   _self.neighborsNodes[v1] = list(set(_self.neighborsNodes[v1]))
   _self.neighborsNodes[v2] = list(set(_self.neighborsNodes[v2]))
   _self.neighborsNodes[v3] = list(set(_self.neighborsNodes[v3]))
   _self.neighborsNodes[v4] = list(set(_self.neighborsNodes[v4]))
   _self.neighborsNodes[v5] = list(set(_self.neighborsNodes[v5]))
   _self.neighborsNodes[v6] = list(set(_self.neighborsNodes[v6]))
   _self.neighborsNodes[v7] = list(set(_self.neighborsNodes[v7]))
   _self.neighborsNodes[v8] = list(set(_self.neighborsNodes[v8]))
   _self.neighborsNodes[v9] = list(set(_self.neighborsNodes[v9]))
   _self.neighborsNodes[v10] = list(set(_self.neighborsNodes[v10]))
   
   _self.neighborsElements[v1].append(e)  
   _self.neighborsElements[v2].append(e)  
   _self.neighborsElements[v3].append(e)  
   _self.neighborsElements[v4].append(e)  
   _self.neighborsElements[v5].append(e)  
   _self.neighborsElements[v6].append(e)  
   _self.neighborsElements[v7].append(e)  
   _self.neighborsElements[v8].append(e)  
   _self.neighborsElements[v9].append(e)  
   _self.neighborsElements[v10].append(e)  
   
   _self.nodes_linear.append(v1)  
   _self.nodes_linear.append(v2)  
   _self.nodes_linear.append(v3)  

   _self.nodes_quad.append(v4)  
   _self.nodes_quad.append(v5)  
   _self.nodes_quad.append(v6)  
   _self.nodes_quad.append(v7)  
   _self.nodes_quad.append(v8)  
   _self.nodes_quad.append(v9)  
   _self.nodes_quad.append(v10)  

   x_a = _self.x[v1] - _self.x[v2]
   x_b = _self.x[v2] - _self.x[v3]
   x_c = _self.x[v3] - _self.x[v1]
   
   y_a = _self.y[v1] - _self.y[v2]
   y_b = _self.y[v2] - _self.y[v3]
   y_c = _self.y[v3] - _self.y[v1]
   
   length1 = np.sqrt(x_a**2 + y_a**2)
   length2 = np.sqrt(x_b**2 + y_b**2)
   length3 = np.sqrt(x_c**2 + y_c**2)

   length.append(length1)
   length.append(length2)
   length.append(length3)
   
  _self.minLengthMesh = min(length)

  _self.nodes_linear = list(set(_self.nodes_linear))  
  _self.nodes_quad = list(set(_self.nodes_quad))


  for i in range(0, _self.numNodes):
   for j in _self.neighborsNodes[i]:
    _self.far_neighborsNodes[i].extend(_self.neighborsNodes[j]) 
    _self.far_neighborsElements[i].extend(_self.neighborsElements[j]) 
 
   _self.far_neighborsNodes[i] = list(set(_self.far_neighborsNodes[i])\
                                     - set(_self.neighborsNodes[i]))
   
   _self.far_neighborsElements[i] = list(set(_self.far_neighborsElements[i])\
                                        - set(_self.neighborsElements[i]))

