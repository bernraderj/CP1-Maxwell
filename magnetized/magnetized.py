from __future__ import print_function
from dolfin import *
from ufl import nabla_div

def main():
    solve('horseshoe_magnet.xml', "horseshoe_magnet_physical_region.xml", "horseshoe_magnet_facet_region.xml", "magnetized_solution_horseshoe_magnet.pvd")
    solve('sphere_magnet.xml', "sphere_magnet_physical_region.xml", "sphere_magnet_facet_region.xml", "sphere_solution_horseshoe_magnet.pvd")
    solve('cube_magnet.xml', "cube_magnet_physical_region.xml", "cube_magnet_facet_region.xml", "magnetized_solution_cube_magnet.pvd")
    

def solve(mesh_file, cell_regions_file, facet_regions_file, output_file):

    # Import mesh with volume and facet regions
    mesh = Mesh(mesh_file)
    cell_regions = MeshFunction("size_t", mesh, cell_regions_file)
    facet_regions = MeshFunction("size_t", mesh, facet_regions_file)

    dx = Measure('dx', domain=mesh, subdomain_data=cell_regions)

    #Define function space
    V = FunctionSpace(mesh, 'P', 1)
    W = VectorFunctionSpace(mesh, 'P', 1)

    #Magnetization
    M0 = Expression('1.0', '0.0', '0.0', degree=0)
    M1 = Expression('0.0', '0.0', '0.0', degree=0)

    # Define variational problem
    u = TrialFunction(V) # scalar potential
    v = TestFunction(V)

    #f0 = div(M0)
    #f1 = div(M1)
    #f0 = M0(0)
    #f0 = Expression('div(M0)', degree=0)
    #f0 = div(M0)
    #f1 = M1(0)
    #TODO: Bugfixing
    f0 = project(div(M0), V)
    f1 = project(div(M1), V)


    u = TrialFunction(V)
    v = TestFunction(V)
    f = Constant(0)
    a = (dot(grad(u), grad(v)))*dx
    #L = f*v*dx(0) + f*v*dx(1)
    L = f0*v*dx(0) + f1*v*dx(1)

    # Dirichlet boundary (1 at (0,0,0))
    bc1 = DirichletBC(V, Constant(1.0), facet_regions, 3)

    #TODO: smaller tol without No-Facets-Warning
    def origin(x, on_boundary):
        #tol = 1E-14
        tol = 5
        return (near(x[0],0,tol) and near(x[1],0,tol) and near(x[2],0,tol))
        #tol = 2E-1
        #return near(x[0], 0, tol) #and near(x[1], 0, tol) and near(x[2], 0, tol)
        #return (near(x[0],0,tol) and near(x[1],0,tol) and near(x[2],0,tol))

        #return near(x[0],0,DOLFIN_EPS)
        #return near(x[0],0,DOLFIN_EPS) and near(x[1],0,DOLFIN_EPS) and near(x[2],0,DOLFIN_EPS)

    bc2 = DirichletBC(V, Constant(1.0), origin)

    # Compute solution
    u = Function(V)
    solve(a == L, u, [bc1,bc2])


    #stray field
    H = project(-grad(u), W)

    #Export
    vtkfile = File(output_file)
    vtkfile << H

main()