import import_vtk

npoints, nelem, IEN, x, y, vx, vy, scalar1, scalar2, scalar3 = import_vtk.vtkfile_linear("marquesleandro/axiHalfPoiseuille/results/psi_boundary_inflow_refinedmesh/psi_boundary_inflow_refinedmesh599.vtk")

print npoints
