from __future__ import print_function
from dolfin import *
from ufl import nabla_div

def main():
    #solver("horseshoe_magnet.xml", "horseshoe_magnet_physical_region.xml", "horseshoe_magnet_facet_region.xml", "magnetized_solution_horseshoe_magnet.pvd", "potential_horseshoe_magnet.pvd")
    #solver("sphere_magnet.xml", "sphere_magnet_physical_region.xml", "sphere_magnet_facet_region.xml", "magnetized_solution_sphere_magnet.pvd", "potential_sphere_magnet.pvd")
    #solver("cube_magnet.xml", "cube_magnet_physical_region.xml", "cube_magnet_facet_region.xml", "magnetized_solution_cube_magnet.pvd", "potential_cube_magnet.pvd")
    solver("cube_magnet_30.xml", "cube_magnet_30_physical_region.xml", "cube_magnet_30_facet_region.xml", "30_magnetized_solution_cube_magnet.pvd", "30_potential_cube_magnet.pvd")
    solver("cube_magnet_50.xml", "cube_magnet_50_physical_region.xml", "cube_magnet_50_facet_region.xml", "50_magnetized_solution_cube_magnet.pvd", "50_potential_cube_magnet.pvd")
    #solver("cube_magnet_70.xml", "cube_magnet_70_physical_region.xml", "cube_magnet_70_facet_region.xml", "70_magnetized_solution_cube_magnet.pvd", "70_potential_cube_magnet.pvd")
    

def solver(mesh_file, cell_regions_file, facet_regions_file, output_file, output_file_potential):

    # Import mesh with volume and facet regions
    mesh = Mesh(mesh_file)
    cell_regions = MeshFunction("size_t", mesh, cell_regions_file)
    facet_regions = MeshFunction("size_t", mesh, facet_regions_file)

    dx = Measure('dx', domain=mesh, subdomain_data=cell_regions)

    #Define function space
    V = FunctionSpace(mesh, 'P', 1)
    W = VectorFunctionSpace(mesh, 'P', 1)

    #Magnetization
    M0 = Expression(('1.0', '0.0', '0.0'), degree=0)
    M1 = Expression(('0.0', '0.0', '0.0'), degree=0)

    # Define variational problem
    u = TrialFunction(V) # scalar potential
    v = TestFunction(V)

    f0 = M0
    f1 = M1

    u = TrialFunction(V)
    v = TestFunction(V)
    a = (dot(grad(u), grad(v)))*dx
    L = (dot(f0, grad(v)))*dx(1) + (dot(f1, grad(v)))*dx(2)

    # Dirichlet boundary (0 at outer boundary)
    bc1 = DirichletBC(V, Constant(0.0), facet_regions, 3)

    # Compute solution
    u = Function(V)
    solve(a == L, u, bc1)
    #solve(a == L, u, bc1, solver_parameters={'linear_solver':'mumps'})

    #stray field
    #H = project(-grad(u), W)
    H = project(-grad(u), W, solver_type="mumps")

    #Export
    vtkfile = File(output_file)
    vtkfile << H

    potentialfile = File(output_file_potential)
    potentialfile << u

main()