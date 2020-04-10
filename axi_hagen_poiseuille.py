# =======================
# Importing the libraries
# =======================

import os
initial_path = os.getcwd()

import sys
directory = './lib_class'
sys.path.insert(0, directory)

from tqdm import tqdm
from time import time

import numpy as np
import scipy.sparse as sps
import scipy.sparse.linalg

import search_file
import import_msh
import assembly
import benchmark_problems
import import_vtk
import semi_lagrangian
import export_vtk
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
directory_save = raw_input(" Enter folder name to save simulations: ")
print' ----------------------------------------------------------------------------\n'



print '\n ------------'
print ' IMPORT MESH:'
print ' ------------'

start_time = time()

# Linear Element
if polynomial_option == 1:
 #mesh_name = 'malha_hagen_poiseuille.msh'
 mesh_name = 'malha_half_poiseuille.msh'
 #mesh_name = 'malha_half_poiseuille20.msh'
 equation_number = 3

 directory = search_file.Find(mesh_name)
 if directory == 'File not found':
  sys.exit()

 msh = import_msh.Linear2D(directory, mesh_name, equation_number)
 msh.coord()
 msh.ien()


# Mini Element
elif polynomial_option == 2:
 mesh_name = 'malha_half_poiseuille.msh'
 equation_number = 3

 directory = search_file.Find(mesh_name)
 if directory == 'File not found':
  sys.exit()

 msh = import_msh.Mini2D(directory, mesh_name, equation_number)
 msh.coord()
 msh.ien()

# Quad Element
elif polynomial_option == 3:
 #mesh_name = 'malha_half_poiseuille_quad.msh'
 mesh_name = 'malha_half_poiseuille_quad1.msh'
 equation_number = 3
 
 directory = search_file.Find(mesh_name)
 if directory == 'File not found':
  sys.exit()

 msh = import_msh.Quad2D(directory, mesh_name, equation_number)
 msh.coord()
 msh.ien()

# Cubic Element
elif polynomial_option == 4:
 mesh_name = 'malha_half_poiseuille_cubic.msh'
 equation_number = 3

 directory = search_file.Find(mesh_name)
 if directory == 'File not found':
  sys.exit()

 msh = import_msh.Cubic2D(directory, mesh_name, equation_number)
 msh.coord()
 msh.ien()



npoints                = msh.npoints
nelem                  = msh.nelem
z                      = msh.x
r                      = msh.y
IEN                    = msh.IEN
neumann_edges          = msh.neumann_edges
dirichlet_pts          = msh.dirichlet_pts
neighbors_nodes        = msh.neighbors_nodes
neighbors_elements     = msh.neighbors_elements
far_neighbors_nodes    = msh.far_neighbors_nodes
far_neighbors_elements = msh.far_neighbors_elements
length_min             = msh.length_min
GL                     = msh.GL
nphysical              = msh.nphysical 


CFL = 0.5
#dt = float(CFL*length_min)
#dt = 0.05
#dt = 0.1
dt = 0.01
Re = 100.0
Sc = 1.0

end_time = time()
import_mesh_time = end_time - start_time
print ' time duration: %.1f seconds \n' %import_mesh_time



print ' ---------'
print ' ASSEMBLY:'
print ' ---------'

start_time = time()

Kzzr, Krrr, Mr, M1r, Mr2, Gr, Gz, Grr, Gzr, polynomial_order = assembly.AxiElement2D(polynomial_option, GL, npoints, nelem, IEN, z, r, gausspoints)


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
 zvelocity_LHS0 = sps.lil_matrix.copy(Mr)
 condition_zvelocity = benchmark_problems.axiHagenPoiseuille(nphysical,npoints,z,r)
 condition_zvelocity.zVelocityProfile_condition(dirichlet_pts[1],zvelocity_LHS0,neighbors_nodes)
 vorticity_ibc = condition_zvelocity.ibc
 benchmark_problem = condition_zvelocity.benchmark_problem

 # Applying vr condition
 rvelocity_LHS0 = sps.lil_matrix.copy(Mr)
 condition_rvelocity = benchmark_problems.axiHagenPoiseuille(nphysical,npoints,z,r)
 condition_rvelocity.rVelocityProfile_condition(dirichlet_pts[2],zvelocity_LHS0,neighbors_nodes)
 
 # Applying psi condition
 streamfunction_LHS0 = sps.lil_matrix.copy(Kzzr) + sps.lil_matrix.copy(Krrr) + 2.0*sps.lil_matrix.copy(Gr)
 condition_streamfunction = benchmark_problems.axiHagenPoiseuille(nphysical,npoints,z,r)
 condition_streamfunction.streamfunction_condition(dirichlet_pts[3],streamfunction_LHS0,neighbors_nodes)

 # Applying vorticity condition
 #condition_vorticity = benchmark_problems.axiHagen_Poiseuille(nphysical,npoints,z,r)
 #condition_vorticity.vorticity_condition(dirichlet_pts[4])
 #vorticity_ibc = condition_vorticity.ibc



# Mini Element
elif polynomial_option == 2:

 # Applying vz condition
 zvelocity_LHS0 = sps.lil_matrix.copy(Mr)
 condition_zvelocity = benchmark_problems.axiHagenPoiseuille(nphysical,npoints,z,r)
 condition_zvelocity.zVelocityProfile_condition(dirichlet_pts[1],zvelocity_LHS0,neighbors_nodes)
 vorticity_ibc = condition_zvelocity.ibc
 benchmark_problem = condition_zvelocity.benchmark_problem

 # Applying vr condition
 rvelocity_LHS0 = sps.lil_matrix.copy(Mr)
 condition_rvelocity = benchmark_problems.axiHagenPoiseuille(nphysical,npoints,z,r)
 condition_rvelocity.rVelocityProfile_condition(dirichlet_pts[2],zvelocity_LHS0,neighbors_nodes)
 
 # Applying psi condition
 streamfunction_LHS0 = sps.lil_matrix.copy(Kzzr) + sps.lil_matrix.copy(Krrr) + 2.0*sps.lil_matrix.copy(Gr)
 condition_streamfunction = benchmark_problems.axiHagenPoiseuille(nphysical,npoints,z,r)
 condition_streamfunction.streamfunction_condition(dirichlet_pts[3],streamfunction_LHS0,neighbors_nodes)

 # Applying vorticity condition
 #condition_vorticity = benchmark_problems.axiHagen_Poiseuille(nphysical,npoints,z,r)
 #condition_vorticity.vorticity_condition(dirichlet_pts[4])
 #vorticity_ibc = condition_vorticity.ibc




# Quad Element
elif polynomial_option == 3:

 # Applying vz condition
 zvelocity_LHS0 = sps.lil_matrix.copy(Mr)
 condition_zvelocity = benchmark_problems.axiQuadHagenPoiseuille(nphysical,npoints,z,r)
 condition_zvelocity.zVelocityProfile_condition(dirichlet_pts[1],zvelocity_LHS0,neighbors_nodes)
 vorticity_ibc = condition_zvelocity.ibc
 benchmark_problem = condition_zvelocity.benchmark_problem

 # Applying vr condition
 rvelocity_LHS0 = sps.lil_matrix.copy(Mr)
 condition_rvelocity = benchmark_problems.axiQuadHagenPoiseuille(nphysical,npoints,z,r)
 condition_rvelocity.rVelocityProfile_condition(dirichlet_pts[2],zvelocity_LHS0,neighbors_nodes)
 
 # Applying psi condition
 streamfunction_LHS0 = sps.lil_matrix.copy(Kzzr) + sps.lil_matrix.copy(Krrr) + 2.0*sps.lil_matrix.copy(Gr)
 condition_streamfunction = benchmark_problems.axiQuadHagenPoiseuille(nphysical,npoints,z,r)
 condition_streamfunction.streamfunction_condition(dirichlet_pts[3],streamfunction_LHS0,neighbors_nodes)

 # Applying vorticity condition
 #condition_vorticity = benchmark_problems.axiQuadHagen_Poiseuille(nphysical,npoints,z,r)
 #condition_vorticity.vorticity_condition(dirichlet_pts[4])
 #vorticity_ibc = condition_vorticity.ibc
# ---------------------------------------------------------------------------------



# -------------------------- Initial condition ------------------------------------
vz = np.copy(condition_zvelocity.bc_1)
vr = np.copy(condition_rvelocity.bc_1)
psi = np.copy(condition_streamfunction.bc_1)
w = np.zeros([npoints,1], dtype = float)
# ---------------------------------------------------------------------------------




#---------- Step 1 - Compute the vorticity and stream field --------------------
# -----Vorticity initial-----
vorticity_RHS = sps.lil_matrix.dot(Gzr,vr) - sps.lil_matrix.dot(Grr,vz)
vorticity_LHS = sps.lil_matrix.copy(Mr)
w = scipy.sparse.linalg.cg(vorticity_LHS,vorticity_RHS,w, maxiter=1.0e+05, tol=1.0e-05)
w = w[0].reshape((len(w[0]),1))


# -----Streamline initial-----
streamfunction_RHS = sps.lil_matrix.dot(Mr2,w)
streamfunction_RHS = np.multiply(streamfunction_RHS,condition_streamfunction.bc_2)
streamfunction_RHS = streamfunction_RHS + condition_streamfunction.bc_dirichlet
psi = scipy.sparse.linalg.cg(condition_streamfunction.LHS,streamfunction_RHS,psi, maxiter=1.0e+05, tol=1.0e-05)
psi = psi[0].reshape((len(psi[0]),1))
#----------------------------------------------------------------------------------




# -------------------------- Import VTK File ------------------------------------
#npoints, nelem, IEN, z, r, vz, vr, w, w, psi = import_vtk.vtkfile_linear("/home/marquesleandro/axiHagenPoiseuille/results/linear8/linear8290.vtk")
#----------------------------------------------------------------------------------




end_time = time()
bc_apply_time = end_time - start_time
print ' time duration: %.1f seconds \n' %bc_apply_time





print ' -----------------------------'
print ' PARAMETERS OF THE SIMULATION:'
print ' -----------------------------'

print ' Mesh: %s' %mesh_name
print ' Number of equation: %s' %equation_number
print ' Number of nodes: %s' %npoints
print ' Number of elements: %s' %nelem
print ' Smallest edge length: %f' %length_min
print ' Time step: %s' %dt
print ' Number of time iteration: %s' %nt
print ' Reynolds number: %s' %Re
print ' Schmidt number: %s' %Sc
print ""


print ' ----------------------------'
print ' SOLVE THE LINEARS EQUATIONS:'
print ' ---------------------------- \n'

print ' Saving simulation in %s \n' %directory_save



start_time = time()
os.chdir(initial_path)



# ------------------------ Export VTK File ---------------------------------------
# Linear and Mini Elements
if polynomial_option == 1 or polynomial_option == 2:   
 save = export_vtk.Linear2D(z,r,IEN,npoints,nelem,w,w,psi,vz,vr)
 save.create_dir(directory_save)
 save.saveVTK(directory_save + str(0))

# Quad Element
elif polynomial_option == 3:   
 save = export_vtk.Quad2D(z,r,IEN,npoints,nelem,w,w,psi,vz,vr)
 save.create_dir(directory_save)
 save.saveVTK(directory_save + str(0))
# ---------------------------------------------------------------------------------

vorticity_bc_1 = np.zeros([npoints,1], dtype = float) 
vz_old = np.zeros([npoints,1], dtype = float)
vr_old = np.zeros([npoints,1], dtype = float)
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
 
 print ' Mesh: %s' %mesh_name
 print ' Number of equation: %s' %equation_number
 print ' Number of nodes: %s' %npoints
 print ' Number of elements: %s' %nelem
 print ' Smallest edge length: %f' %length_min
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

 Mv = assembly.AxiAssembleMv(polynomial_option, GL, npoints, nelem, IEN, z, r, vr, gausspoints)
 print ""
 # ---------------------------------------------------------------------------------



 # ------------------------ SOLVE LINEAR EQUATIONS ----------------------------------
 print ' ----------------------------'
 print ' SOLVE THE LINEARS EQUATIONS:'
 print ' ----------------------------'
 print ""
 print ' Saving simulation in %s' %directory_save
 print ""

 start_solver_time = time()


 #---------- Step 2 - Compute the boundary conditions for vorticity --------------
 vorticity_RHS = sps.lil_matrix.dot(Gzr,vr) - sps.lil_matrix.dot(Grr,vz)
 vorticity_LHS = sps.lil_matrix.copy(Mr)
 vorticity_bc_1 = scipy.sparse.linalg.cg(vorticity_LHS,vorticity_RHS,vorticity_bc_1, maxiter=1.0e+05, tol=1.0e-05)
 vorticity_bc_1 = vorticity_bc_1[0].reshape((len(vorticity_bc_1[0]),1))

 if polynomial_option == 1 or polynomial_option == 2:   
  for i in range(0,len(dirichlet_pts[2])):
   line = dirichlet_pts[2][i][0]
   v1 = dirichlet_pts[2][i][1] - 1
   v2 = dirichlet_pts[2][i][2] - 1

   if line == 8:
    vorticity_bc_1[v1] = 0.0
    vorticity_bc_1[v2] = 0.0

    np.append(vorticity_ibc,v1)
    np.append(vorticity_ibc,v2)

  vorticity_ibc = np.unique(vorticity_ibc)

 elif polynomial_option == 3:
  for i in range(0,len(dirichlet_pts[2])):
   line = dirichlet_pts[2][i][0]
   v1 = dirichlet_pts[2][i][1] - 1
   v2 = dirichlet_pts[2][i][2] - 1
   v3 = dirichlet_pts[2][i][3] - 1

   if line == 8:
    vorticity_bc_1[v1] = 0.0
    vorticity_bc_1[v2] = 0.0
    vorticity_bc_1[v3] = 0.0

    np.append(vorticity_ibc,v1)
    np.append(vorticity_ibc,v2)
    np.append(vorticity_ibc,v3)

  vorticity_ibc = np.unique(vorticity_ibc)

 # Gaussian elimination
 vorticity_bc_dirichlet = np.zeros([npoints,1], dtype = float)
 vorticity_bc_neumann = np.zeros([npoints,1], dtype = float)
 vorticity_bc_2 = np.ones([npoints,1], dtype = float)

 vorticity_LHS = (np.copy(Mr)/dt) + (1.0/Re)*np.copy(Krrr) + (1.0/Re)*np.copy(Kzzr)  + (1.0/Re)*np.copy(M1r) - Mv
 for mm in vorticity_ibc:
  for nn in neighbors_nodes[mm]:
   vorticity_bc_dirichlet[nn] -= float(vorticity_LHS[nn,mm]*vorticity_bc_1[mm])
   vorticity_LHS[nn,mm] = 0.0
   vorticity_LHS[mm,nn] = 0.0
   
  vorticity_LHS[mm,mm] = 1.0
  vorticity_bc_dirichlet[mm] = vorticity_bc_1[mm]
  vorticity_bc_2[mm] = 0.0
 #----------------------------------------------------------------------------------



 #---------- Step 3 - Solve the vorticity transport equation ----------------------
 # Taylor Galerkin Scheme
 if scheme_option == 1:
  scheme_name = 'Taylor Galerkin'
  A = np.copy(M)/dt 
  vorticity_RHS = sps.lil_matrix.dot(A,w) - np.multiply(vz,sps.lil_matrix.dot(Gz,w))\
        - np.multiply(vr,sps.lil_matrix.dot(Gr,w))\
        - (dt/2.0)*np.multiply(vz,(np.multiply(vz,sps.lil_matrix.dot(Kzz,w)) + np.multiply(vr,sps.lil_matrix.dot(Krz,w))))\
        - (dt/2.0)*np.multiply(vr,(np.multiply(vz,sps.lil_matrix.dot(Kzr,w)) + np.multiply(vr,sps.lil_matrix.dot(Krr,w))))
  vorticity_RHS = np.multiply(vorticity_RHS,vorticity_bc_2)
  vorticity_RHS = vorticity_RHS + vorticity_bc_dirichlet
  w = scipy.sparse.linalg.cg(vorticity_LHS,vorticity_RHS,w, maxiter=1.0e+05, tol=1.0e-05)
  w = w[0].reshape((len(w[0]),1))



 # Semi-Lagrangian Scheme
 elif scheme_option == 2:

  # Linear Element   
  if polynomial_option == 1:
   scheme_name = 'Semi Lagrangian Linear'

   w_d = semi_lagrangian.Linear2D(npoints, neighbors_elements, IEN, z, r, vz, vr, dt, w)

   A = np.copy(Mr)/dt
   vorticity_RHS = sps.lil_matrix.dot(A,w_d)

   vorticity_RHS = vorticity_RHS + (1.0/Re)*vorticity_bc_neumann
   vorticity_RHS = np.multiply(vorticity_RHS,vorticity_bc_2)
   vorticity_RHS = vorticity_RHS + vorticity_bc_dirichlet

   w = scipy.sparse.linalg.cg(vorticity_LHS,vorticity_RHS,w, maxiter=1.0e+05, tol=1.0e-05)
   w = w[0].reshape((len(w[0]),1))



  # Mini Element   
  elif polynomial_option == 2:
   scheme_name = 'Semi Lagrangian Mini'

   w_d = semi_lagrangian.Mini2D(npoints, neighbors_elements, IEN, z, r, vz, vr, dt, w)

   A = np.copy(Mr)/dt
   vorticity_RHS = sps.lil_matrix.dot(A,w_d)

   vorticity_RHS = vorticity_RHS + (1.0/Re)*vorticity_bc_neumann
   vorticity_RHS = np.multiply(vorticity_RHS,vorticity_bc_2)
   vorticity_RHS = vorticity_RHS + vorticity_bc_dirichlet

   w = scipy.sparse.linalg.cg(vorticity_LHS,vorticity_RHS,w, maxiter=1.0e+05, tol=1.0e-05)
   w = w[0].reshape((len(w[0]),1))



  # Quad Element   
  elif polynomial_option == 3:
   scheme_name = 'Semi Lagrangian Quad'

   w_d = semi_lagrangian.Quad2D(npoints, neighbors_elements, IEN, z, r, vz, vr, dt, w)

   A = np.copy(Mr)/dt
   vorticity_RHS = sps.lil_matrix.dot(A,w_d)

   vorticity_RHS = vorticity_RHS + (1.0/Re)*vorticity_bc_neumann
   vorticity_RHS = np.multiply(vorticity_RHS,vorticity_bc_2)
   vorticity_RHS = vorticity_RHS + vorticity_bc_dirichlet

   w = scipy.sparse.linalg.cg(vorticity_LHS,vorticity_RHS,w, maxiter=1.0e+05, tol=1.0e-05)
   w = w[0].reshape((len(w[0]),1)) 
 #----------------------------------------------------------------------------------



 #---------- Step 4 - Solve the streamline equation --------------------------------
 # Solve Streamline
 # psi condition
 streamfunction_RHS = sps.lil_matrix.dot(Mr2,w)
 streamfunction_RHS = np.multiply(streamfunction_RHS,condition_streamfunction.bc_2)
 streamfunction_RHS = streamfunction_RHS + condition_streamfunction.bc_dirichlet
 psi = scipy.sparse.linalg.cg(condition_streamfunction.LHS,streamfunction_RHS,psi, maxiter=1.0e+05, tol=1.0e-05)
 psi = psi[0].reshape((len(psi[0]),1))
 #----------------------------------------------------------------------------------



 #---------- Step 5 - Compute the velocity field -----------------------------------
 # Velocity vz
 vz_old = np.copy(vz)
 zvelocity_RHS = sps.lil_matrix.dot(Gr,psi)
 zvelocity_RHS = np.multiply(zvelocity_RHS,condition_zvelocity.bc_2)
 zvelocity_RHS = zvelocity_RHS + condition_zvelocity.bc_dirichlet
 vz = scipy.sparse.linalg.cg(condition_zvelocity.LHS,zvelocity_RHS,vz, maxiter=1.0e+05, tol=1.0e-05)
 vz = vz[0].reshape((len(vz[0]),1))
 
 # Velocity vr
 vr_old = np.copy(vr)
 rvelocity_RHS = -sps.lil_matrix.dot(Gz,psi)
 rvelocity_RHS = np.multiply(rvelocity_RHS,condition_rvelocity.bc_2)
 rvelocity_RHS = rvelocity_RHS + condition_rvelocity.bc_dirichlet
 vr = scipy.sparse.linalg.cg(condition_rvelocity.LHS,rvelocity_RHS,vr, maxiter=1.0e+05, tol=1.0e-05)
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
  save = export_vtk.Linear2D(z,r,IEN,npoints,nelem,w,w,psi,vz,vr)
  save.create_dir(directory_save)
  save.saveVTK(directory_save + str(t))

 # Quad Element
 elif polynomial_option == 3:   
  save = export_vtk.Quad2D(z,r,IEN,npoints,nelem,w,w,psi,vz,vr)
  save.create_dir(directory_save)
  save.saveVTK(directory_save + str(t))


 end_time = time()
 export_time_solver = end_time - start_time
 print ' time duration: %.1f seconds' %export_time_solver
 print ""
 # ---------------------------------------------------------------------------------




 # ------------------------ CHECK STEADY STATE ----------------------------------
 #vz_dif = np.sqrt((vz-vz_old)**2)
 #vr_dif = np.sqrt((vr-vr_old)**2)
 #if np.all(vz_dif < 1e-50) and np.all(vr_dif < 1e-50):
 # end_type = 1
 # break
 # ---------------------------------------------------------------------------------

 # ------------------------ CHECK CONVERGENCE RESULT ----------------------------------
 if np.linalg.norm(vz) > 10e2 or np.linalg.norm(vr) > 10e2:
  end_type = 2
  break
 # ---------------------------------------------------------------------------------
 






end_time = time()
solution_time = end_time - start_time
print ' time duration: %.1f seconds \n' %solution_time


print ' ----------------'
print ' SAVING RELATORY:'
print ' ----------------'
print ""

if end_type == 0:
 print ' END SIMULATION. NOT STEADY STATE'
 print ' Relatory saved in %s' %directory_save
 print ""

elif end_type == 1:
 print ' END SIMULATION. STEADY STATE'
 print ' Relatory saved in %s' %directory_save
 print ""

elif end_type == 2:
 print ' END SIMULATION. ERROR CONVERGENCE RESULT'
 print ' Relatory saved in %s' %directory_save
 print ""




# -------------------------------- Export Relatory ---------------------------------------
relatory.export(save.path, directory_save, sys.argv[0], benchmark_problem, scheme_name, mesh_name, equation_number, npoints, nelem, length_min, dt, nt, Re, Sc, import_mesh_time, assembly_time, bc_apply_time, solution_time, polynomial_order, gausspoints)
# ----------------------------------------------------------------------------------------



