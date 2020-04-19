# =======================
# Importing the libraries
# =======================

import os
initial_path = os.getcwd()

import sys
folderClass = './libClass'
sys.path.insert(0, folderClass)

from tqdm import tqdm
from time import time

import numpy as np
import scipy.sparse as sps
import scipy.sparse.linalg

import searchMSH
import importMSH
import assembly
import benchmarkProblems
import importVTK
import semiLagrangian
import exportVTK
import relatory



print '''
               COPYRIGHT                    
 ======================================
 Simulator: %s
 created by Leandro Marques at 02/2019
 e-mail: marquesleandro67@gmail.com
 Gesar Search Group
 State University of the Rio de Janeiro
 ======================================
\n''' %sys.argv[0]



print ' ------'
print ' INPUT:'
print ' ------'

print ""
print ' ----------------------------------------------------------------------------'
print ' (1) - Linear Element'
print ' (2) - Mini Element'
print ' (3) - Quadratic Element'
print ' (4) - Cubic Element'
polynomial_option = int(raw_input("\n Enter polynomial degree option above: "))
print' ----------------------------------------------------------------------------\n'


print ' ----------------------------------------------------------------------------'
print ' 3 Gauss Points'
print ' 4 Gauss Points'
print ' 6 Gauss Points'
print ' 12 Gauss Points'
gausspoints = int(raw_input("\n Enter Gauss Points Number option above: "))
print' ----------------------------------------------------------------------------\n'


print ' ----------------------------------------------------------------------------'
print ' (1) - Taylor Galerkin Scheme'
print ' (2) - Semi Lagrangian Scheme'
scheme_option = int(raw_input("\n Enter simulation scheme option above: "))
print' ----------------------------------------------------------------------------\n'


print ' ----------------------------------------------------------------------------'
nt = int(raw_input(" Enter number of time interations (nt): "))
print' ----------------------------------------------------------------------------\n'


print ' ----------------------------------------------------------------------------'
folderResults = raw_input(" Enter folder name to save simulations: ")
print' ----------------------------------------------------------------------------\n'



print '\n ------------'
print ' IMPORT MESH:'
print ' ------------'

start_time = time()

# Linear Element
if polynomial_option == 1:
 mshFileName = 'linearHagenPoiseuille.msh'

 pathMSHFile = searchMSH.Find(mshFileName)
 if pathMSHFile == 'File not found':
  sys.exit()

 mesh = importMSH.Linear2D(pathMSHFile, mshFileName)


# Mini Element
elif polynomial_option == 2:
 mshFileName = 'linearHagenPoiseuille.msh'
 equation_number = 3

 pathMSHFile = searchMSH.Find(mshFileName)
 if pathMSHFile == 'File not found':
  sys.exit()

 mesh = importMSH.Mini2D(pathMSHFile, mshFileName, equation_number)
 mesh.coord()
 mesh.ien()

# Quad Element
elif polynomial_option == 3:
 mshFileName = 'quadHagenPoiseuille.msh'
 equation_number = 3
 
 pathMSHFile = searchMSH.Find(mshFileName)
 if pathMSHFile == 'File not found':
  sys.exit()

 mesh = importMSH.Quad2D(pathMSHFile, mshFileName, equation_number)
 mesh.coord()
 mesh.ien()

# Cubic Element
elif polynomial_option == 4:
 mshFileName = 'cubicHagenPoiseuille.msh'
 equation_number = 3

 pathMSHFile = searchMSH.Find(mshFileName)
 if pathMSHFile == 'File not found':
  sys.exit()

 mesh = importMSH.Cubic2D(pathMSHFile, mshFileName, equation_number)
 mesh.coord()
 mesh.ien()



numNodes               = mesh.numNodes
numElements            = mesh.numElements
z                      = mesh.x
r                      = mesh.y
IEN                    = mesh.IEN
boundaryEdges          = mesh.boundaryEdges
boundaryNodes          = mesh.boundaryNodes
neighborsNodes         = mesh.neighborsNodes
neighborsElements      = mesh.neighborsElements
minLengthMesh          = mesh.minLengthMesh
FreedomDegree          = mesh.FreedomDegree
numPhysical            = mesh.numPhysical 


CFL = 0.5
#dt = float(CFL*minLengthMesh)
dt = 0.1   #linear result ok 
#dt = 0.05  #quad
Re = 100.0
Sc = 1.0

end_time = time()
import_mesh_time = end_time - start_time
print ' time duration: %.1f seconds \n' %import_mesh_time



print ' ---------'
print ' ASSEMBLY:'
print ' ---------'

start_time = time()

Kzzr, Krrr, Mr, M1r, Mr2, Gr, Gz, Grr, Gzr, polynomial_order = assembly.AxiElement2D(polynomial_option, FreedomDegree, numNodes, numElements, IEN, z, r, gausspoints)


end_time = time()
assembly_time = end_time - start_time
print ' time duration: %.1f seconds \n' %assembly_time




print ' --------------------------------'
print ' INITIAL AND BOUNDARY CONDITIONS:'
print ' --------------------------------'

start_time = time()


# ------------------------ Boundaries Conditions ----------------------------------

# Linear Element
if polynomial_option == 1:

 # Applying vz condition
 zVelocityLHS0 = sps.lil_matrix.copy(Mr)
 zVelocityBC = benchmarkProblems.axiHagenPoiseuille(numPhysical,numNodes,z,r)
 zVelocityBC.zVelocityCondition(boundaryEdges,zVelocityLHS0,neighborsNodes)
 benchmark_problem = zVelocityBC.benchmark_problem

 # Applying vr condition
 rVelocityLHS0 = sps.lil_matrix.copy(Mr)
 rVelocityBC = benchmarkProblems.axiHagenPoiseuille(numPhysical,numNodes,z,r)
 rVelocityBC.rVelocityCondition(boundaryEdges,zVelocityLHS0,neighborsNodes)
 
 # Applying psi condition
 streamFunctionLHS0 = sps.lil_matrix.copy(Kzzr) + sps.lil_matrix.copy(Krrr) + 2.0*sps.lil_matrix.copy(Gr)
 streamFunctionBC = benchmarkProblems.axiHagenPoiseuille(numPhysical,numNodes,z,r)
 streamFunctionBC.streamFunctionCondition(boundaryEdges,streamFunctionLHS0,neighborsNodes)

 # Applying vorticity condition
 vorticityDirichletNodes = boundaryNodes



# Mini Element
elif polynomial_option == 2:

 # Applying vz condition
 zVelocityLHS0 = sps.lil_matrix.copy(Mr)
 zVelocityBC = benchmarkProblems.axiHagenPoiseuille(numPhysical,numNodes,z,r)
 zVelocityBC.zVelocityProfile_condition(boundaryEdges[1],zVelocityLHS0,neighborsNodes)
 vorticityDirichletNodes = zVelocityBC.dirichletNodes
 benchmark_problem = zVelocityBC.benchmark_problem

 # Applying vr condition
 rVelocityLHS0 = sps.lil_matrix.copy(Mr)
 rVelocityBC = benchmarkProblems.axiHagenPoiseuille(numPhysical,numNodes,z,r)
 rVelocityBC.rVelocityProfile_condition(boundaryEdges[2],zVelocityLHS0,neighborsNodes)
 
 # Applying psi condition
 streamFunctionLHS0 = sps.lil_matrix.copy(Kzzr) + sps.lil_matrix.copy(Krrr) + 2.0*sps.lil_matrix.copy(Gr)
 streamFunctionBC = benchmarkProblems.axiHagenPoiseuille(numPhysical,numNodes,z,r)
 streamFunctionBC.streamFunctioncondition(boundaryEdges[3],streamFunctionLHS0,neighborsNodes)

 # Applying vorticity condition
 #condition_vorticity = benchmarkProblems.axiHagen_Poiseuille(numPhysical,numNodes,z,r)
 #condition_vorticity.vorticitycondition(boundaryEdges[4])
 #vorticityDirichletNodes = condition_vorticity.dirichletNodes




# Quad Element
elif polynomial_option == 3:

 # Applying vz condition
 zVelocityLHS0 = sps.lil_matrix.copy(Mr)
 zVelocityBC = benchmarkProblems.axiQuadHagenPoiseuille(numPhysical,numNodes,z,r)
 zVelocityBC.zVelocityProfile_condition(boundaryEdges[1],zVelocityLHS0,neighborsNodes)
 vorticityDirichletNodes = zVelocityBC.dirichletNodes
 benchmark_problem = zVelocityBC.benchmark_problem

 # Applying vr condition
 rVelocityLHS0 = sps.lil_matrix.copy(Mr)
 rVelocityBC = benchmarkProblems.axiQuadHagenPoiseuille(numPhysical,numNodes,z,r)
 rVelocityBC.rVelocityProfile_condition(boundaryEdges[2],zVelocityLHS0,neighborsNodes)
 
 # Applying psi condition
 streamFunctionLHS0 = sps.lil_matrix.copy(Kzzr) + sps.lil_matrix.copy(Krrr) + 2.0*sps.lil_matrix.copy(Gr)
 streamFunctionBC = benchmarkProblems.axiQuadHagenPoiseuille(numPhysical,numNodes,z,r)
 streamFunctionBC.streamFunctioncondition(boundaryEdges[3],streamFunctionLHS0,neighborsNodes)

 # Applying vorticity condition
 #condition_vorticity = benchmarkProblems.axiQuadHagen_Poiseuille(numPhysical,numNodes,z,r)
 #condition_vorticity.vorticitycondition(boundaryEdges[4])
 #vorticityDirichletNodes = condition_vorticity.dirichletNodes
# ---------------------------------------------------------------------------------



# -------------------------- Initial condition ------------------------------------
vz = np.copy(zVelocityBC.aux1BC)
vr = np.copy(rVelocityBC.aux1BC)
psi = np.copy(streamFunctionBC.aux1BC)
w = np.zeros([numNodes,1], dtype = float)
# ---------------------------------------------------------------------------------




#---------- Step 1 - Compute the vorticity and stream field --------------------
# -----Vorticity initial-----
vorticityRHS = sps.lil_matrix.dot(Gzr,vr) - sps.lil_matrix.dot(Grr,vz)
vorticityLHS = sps.lil_matrix.copy(Mr)
w = scipy.sparse.linalg.cg(vorticityLHS,vorticityRHS,w, maxiter=1.0e+05, tol=1.0e-05)
w = w[0].reshape((len(w[0]),1))


# -----Streamline initial-----
streamFunctionRHS = sps.lil_matrix.dot(Mr2,w)
streamFunctionRHS = np.multiply(streamFunctionRHS,streamFunctionBC.aux2BC)
streamFunctionRHS = streamFunctionRHS + streamFunctionBC.dirichletVector
psi = scipy.sparse.linalg.cg(streamFunctionBC.LHS,streamFunctionRHS,psi, maxiter=1.0e+05, tol=1.0e-05)
psi = psi[0].reshape((len(psi[0]),1))
#----------------------------------------------------------------------------------




# -------------------------- Import VTK File ------------------------------------
#numNodes, numElements, IEN, z, r, vz, vr, w, w, psi = importVTK.vtkfile_linear("/home/marquesleandro/axiHagenPoiseuille/results/linear8/linear8290.vtk")
#----------------------------------------------------------------------------------




end_time = time()
bc_apply_time = end_time - start_time
print ' time duration: %.1f seconds \n' %bc_apply_time





print ' -----------------------------'
print ' PARAMETERS OF THE SIMULATION:'
print ' -----------------------------'

print ' Benchmark Problem: %s' %benchmark_problem
print ' Scheme: %s' %str(scheme_option)
print ' Element Type: %s' %str(polynomial_order)
print ' Gaussian Quadrature (Gauss Points): %s' %str(gausspoints)
print ' Mesh: %s' %mshFileName
print ' Number of nodes: %s' %numNodes
print ' Number of elements: %s' %numElements
print ' Smallest edge length: %f' %minLengthMesh
print ' Time step: %s' %dt
print ' Number of time iteration: %s' %nt
print ' Reynolds number: %s' %Re
print ' Schmidt number: %s' %Sc
print ""


print ' ----------------------------'
print ' SOLVE THE LINEARS EQUATIONS:'
print ' ---------------------------- \n'

print ' Saving simulation in %s \n' %folderResults



solution_start_time = time()
os.chdir(initial_path)



# ------------------------ Export VTK File ---------------------------------------
# Linear and Mini Elements
if polynomial_option == 1 or polynomial_option == 2:   
 save = exportVTK.Linear2D(z,r,IEN,numNodes,numElements,w,w,psi,vz,vr)
 save.create_dir(folderResults)
 save.saveVTK(folderResults + str(0))

# Quad Element
elif polynomial_option == 3:   
 save = exportVTK.Quad2D(z,r,IEN,numNodes,numElements,w,w,psi,vz,vr)
 save.create_dir(folderResults)
 save.saveVTK(folderResults + str(0))
# ---------------------------------------------------------------------------------

vorticityAux1BC = np.zeros([numNodes,1], dtype = float) 
vz_old = np.zeros([numNodes,1], dtype = float)
vr_old = np.zeros([numNodes,1], dtype = float)
end_type = 0
for t in tqdm(range(1, nt)):
 print ""
 print '''
                COPYRIGHT                    
  ======================================
  Simulator: %s
  created by Leandro Marques at 02/2019
  e-mail: marquesleandro67@gmail.com
  Gesar Search Group
  State University of the Rio de Janeiro
  ======================================
 ''' %sys.argv[0]



 print ' -----------------------------'
 print ' PARAMETERS OF THE SIMULATION:'
 print ' -----------------------------'

 print ' Benchmark Problem: %s' %benchmark_problem
 print ' Scheme: %s' %str(scheme_option)
 print ' Element Type: %s' %str(polynomial_order)
 print ' Gaussian Quadrature (Gauss Points): %s' %str(gausspoints)
 print ' Mesh: %s' %mshFileName
 print ' Number of nodes: %s' %numNodes
 print ' Number of elements: %s' %numElements
 print ' Smallest edge length: %f' %minLengthMesh
 print ' Time step: %s' %dt
 print ' Number of time iteration: %s' %t
 print ' Reynolds number: %s' %Re
 print ' Schmidt number: %s' %Sc
 print ""



 # ------------------------- ASSEMBLY Mv -------------------------------------------
 print ""
 print ' ------------'
 print ' ASSEMBLY Mv:'
 print ' ------------'

 Mv = assembly.AxiAssembleMv(polynomial_option, FreedomDegree, numNodes, numElements, IEN, z, r, vr, gausspoints)
 print ""
 # ---------------------------------------------------------------------------------



 # ------------------------ SOLVE LINEAR EQUATIONS ----------------------------------
 print ' ----------------------------'
 print ' SOLVE THE LINEARS EQUATIONS:'
 print ' ----------------------------'
 print ""
 print ' Saving simulation in %s' %folderResults
 print ""

 start_solver_time = time()


 #---------- Step 2 - Compute the boundary conditions for vorticity --------------
 vorticityRHS = sps.lil_matrix.dot(Gzr,vr) - sps.lil_matrix.dot(Grr,vz)
 vorticityLHS = sps.lil_matrix.copy(Mr)
 vorticityAux1BC = scipy.sparse.linalg.cg(vorticityLHS,vorticityRHS,vorticityAux1BC, maxiter=1.0e+05, tol=1.0e-05)
 vorticityAux1BC = vorticityAux1BC[0].reshape((len(vorticityAux1BC[0]),1))

 # Gaussian elimination
 vorticityDirichletVector = np.zeros([numNodes,1], dtype = float)
 vorticityNeumannVector = np.zeros([numNodes,1], dtype = float)
 vorticityAux2BC = np.ones([numNodes,1], dtype = float)

 vorticityLHS = (np.copy(Mr)/dt) + (1.0/Re)*np.copy(Krrr) + (1.0/Re)*np.copy(Kzzr)  + (1.0/Re)*np.copy(M1r) - Mv
 for mm in vorticityDirichletNodes:
  for nn in neighborsNodes[mm]:
   vorticityDirichletVector[nn] -= float(vorticityLHS[nn,mm]*vorticityAux1BC[mm])
   vorticityLHS[nn,mm] = 0.0
   vorticityLHS[mm,nn] = 0.0
   
  vorticityLHS[mm,mm] = 1.0
  vorticityDirichletVector[mm] = vorticityAux1BC[mm]
  vorticityAux2BC[mm] = 0.0
 #----------------------------------------------------------------------------------



 #---------- Step 3 - Solve the vorticity transport equation ----------------------
 # Taylor Galerkin Scheme
 if scheme_option == 1:
  scheme_name = 'Taylor Galerkin'
  A = np.copy(M)/dt 
  vorticityRHS = sps.lil_matrix.dot(A,w) - np.multiply(vz,sps.lil_matrix.dot(Gz,w))\
        - np.multiply(vr,sps.lil_matrix.dot(Gr,w))\
        - (dt/2.0)*np.multiply(vz,(np.multiply(vz,sps.lil_matrix.dot(Kzz,w)) + np.multiply(vr,sps.lil_matrix.dot(Krz,w))))\
        - (dt/2.0)*np.multiply(vr,(np.multiply(vz,sps.lil_matrix.dot(Kzr,w)) + np.multiply(vr,sps.lil_matrix.dot(Krr,w))))
  vorticityRHS = np.multiply(vorticityRHS,vorticityAux2BC)
  vorticityRHS = vorticityRHS + vorticityDirichletVector
  w = scipy.sparse.linalg.cg(vorticityLHS,vorticityRHS,w, maxiter=1.0e+05, tol=1.0e-05)
  w = w[0].reshape((len(w[0]),1))



 # Semi-Lagrangian Scheme
 elif scheme_option == 2:

  # Linear Element   
  if polynomial_option == 1:
   scheme_name = 'Semi Lagrangian Linear'

   w_d = semiLagrangian.Linear2D(numNodes, neighborsElements, IEN, z, r, vz, vr, dt, w)

   A = np.copy(Mr)/dt
   vorticityRHS = sps.lil_matrix.dot(A,w_d)

   vorticityRHS = vorticityRHS + (1.0/Re)*vorticityNeumannVector
   vorticityRHS = np.multiply(vorticityRHS,vorticityAux2BC)
   vorticityRHS = vorticityRHS + vorticityDirichletVector

   w = scipy.sparse.linalg.cg(vorticityLHS,vorticityRHS,w, maxiter=1.0e+05, tol=1.0e-05)
   w = w[0].reshape((len(w[0]),1))



  # Mini Element   
  elif polynomial_option == 2:
   scheme_name = 'Semi Lagrangian Mini'

   w_d = semiLagrangian.Mini2D(numNodes, neighborsElements, IEN, z, r, vz, vr, dt, w)

   A = np.copy(Mr)/dt
   vorticityRHS = sps.lil_matrix.dot(A,w_d)

   vorticityRHS = vorticityRHS + (1.0/Re)*vorticityNeumannVector
   vorticityRHS = np.multiply(vorticityRHS,vorticityAux2BC)
   vorticityRHS = vorticityRHS + vorticityDirichletVector

   w = scipy.sparse.linalg.cg(vorticityLHS,vorticityRHS,w, maxiter=1.0e+05, tol=1.0e-05)
   w = w[0].reshape((len(w[0]),1))



  # Quad Element   
  elif polynomial_option == 3:
   scheme_name = 'Semi Lagrangian Quad'

   w_d = semiLagrangian.Quad2D(numNodes, neighborsElements, IEN, z, r, vz, vr, dt, w)

   A = np.copy(Mr)/dt
   vorticityRHS = sps.lil_matrix.dot(A,w_d)

   vorticityRHS = vorticityRHS + (1.0/Re)*vorticityNeumannVector
   vorticityRHS = np.multiply(vorticityRHS,vorticityAux2BC)
   vorticityRHS = vorticityRHS + vorticityDirichletVector

   w = scipy.sparse.linalg.cg(vorticityLHS,vorticityRHS,w, maxiter=1.0e+05, tol=1.0e-05)
   w = w[0].reshape((len(w[0]),1)) 
 #----------------------------------------------------------------------------------



 #---------- Step 4 - Solve the streamline equation --------------------------------
 # Solve Streamline
 # psi condition
 streamFunctionRHS = sps.lil_matrix.dot(Mr2,w)
 streamFunctionRHS = np.multiply(streamFunctionRHS,streamFunctionBC.aux2BC)
 streamFunctionRHS = streamFunctionRHS + streamFunctionBC.dirichletVector
 psi = scipy.sparse.linalg.cg(streamFunctionBC.LHS,streamFunctionRHS,psi, maxiter=1.0e+05, tol=1.0e-05)
 psi = psi[0].reshape((len(psi[0]),1))
 #----------------------------------------------------------------------------------



 #---------- Step 5 - Compute the velocity field -----------------------------------
 # Velocity vz
 vz_old = np.copy(vz)
 zVelocityRHS = sps.lil_matrix.dot(Gr,psi)
 zVelocityRHS = np.multiply(zVelocityRHS,zVelocityBC.aux2BC)
 zVelocityRHS = zVelocityRHS + zVelocityBC.dirichletVector
 vz = scipy.sparse.linalg.cg(zVelocityBC.LHS,zVelocityRHS,vz, maxiter=1.0e+05, tol=1.0e-05)
 vz = vz[0].reshape((len(vz[0]),1))
 
 # Velocity vr
 vr_old = np.copy(vr)
 rVelocityRHS = -sps.lil_matrix.dot(Gz,psi)
 rVelocityRHS = np.multiply(rVelocityRHS,rVelocityBC.aux2BC)
 rVelocityRHS = rVelocityRHS + rVelocityBC.dirichletVector
 vr = scipy.sparse.linalg.cg(rVelocityBC.LHS,rVelocityRHS,vr, maxiter=1.0e+05, tol=1.0e-05)
 vr = vr[0].reshape((len(vr[0]),1))
 #----------------------------------------------------------------------------------

 end_solver_time = time()
 solver_time = end_solver_time - start_solver_time
 print ' time duration: %.1f seconds' %solver_time
 print ""
 #----------------------------------------------------------------------------------
 


 # ------------------------ Export VTK File ---------------------------------------
 print ' ----------------'
 print ' EXPORT VTK FILE:'
 print ' ----------------'


 start_time = time()


 # Linear and Mini Elements
 if polynomial_option == 1 or polynomial_option == 2:   
  save = exportVTK.Linear2D(z,r,IEN,numNodes,numElements,w,w,psi,vz,vr)
  save.create_dir(folderResults)
  save.saveVTK(folderResults + str(t))

 # Quad Element
 elif polynomial_option == 3:   
  save = exportVTK.Quad2D(z,r,IEN,numNodes,numElements,w,w,psi,vz,vr)
  save.create_dir(folderResults)
  save.saveVTK(folderResults + str(t))


 end_time = time()
 export_time_solver = end_time - start_time
 print ' time duration: %.1f seconds' %export_time_solver
 print ""
 # ---------------------------------------------------------------------------------




 # ------------------------ CHECK STEADY STATE ----------------------------------
 if np.all(vz == vz_old) and np.all(vr == vr_old):
  end_type = 1
  break
 # ---------------------------------------------------------------------------------

 # ------------------------ CHECK CONVERGENCE RESULT ----------------------------------
 if np.linalg.norm(vz) > 10e2 or np.linalg.norm(vr) > 10e2:
  end_type = 2
  break
 # ---------------------------------------------------------------------------------
 






end_time = time()
solution_time = end_time - solution_start_time
print ' time duration: %.1f seconds \n' %solution_time


print ' ----------------'
print ' SAVING RELATORY:'
print ' ----------------'
print ""

if end_type == 0:
 print ' END SIMULATION. NOT STEADY STATE'
 print ' Relatory saved in %s' %folderResults
 print ""

elif end_type == 1:
 print ' END SIMULATION. STEADY STATE'
 print ' Relatory saved in %s' %folderResults
 print ""

elif end_type == 2:
 print ' END SIMULATION. ERROR CONVERGENCE RESULT'
 print ' Relatory saved in %s' %folderResults
 print ""




# -------------------------------- Export Relatory ---------------------------------------
relatory.export(save.path, folderResults, sys.argv[0], benchmark_problem, scheme_name, mshFileName, numNodes, numElements, minLengthMesh, dt, nt, Re, Sc, import_mesh_time, assembly_time, bc_apply_time, solution_time, polynomial_order, gausspoints)
# ----------------------------------------------------------------------------------------



