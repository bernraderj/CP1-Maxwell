from __future__ import print_function
from dolfin import *

def convert(mesh_file, cell_regions_file, facet_regions_file, output_file):
    mesh = Mesh(mesh_file)
    cell_regions = MeshFunction("size_t", mesh, cell_regions_file)

    vtkfile = File(output_file)
    vtkfile << cell_regions


convert("horseshoe_magnet.xml", "horseshoe_magnet_physical_region.xml", "horseshoe_magnet_facet_region.xml", "domains_horseshoe.pvd")
convert("sphere_magnet.xml", "sphere_magnet_physical_region.xml", "sphere_magnet_facet_region.xml", "domains_sphere.pvd")
convert("cube_magnet.xml", "cube_magnet_physical_region.xml", "cube_magnet_facet_region.xml", "domains_cube.pvd")
