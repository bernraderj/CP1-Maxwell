from fenics import *
import pygmsh
import os

class GenerateMesh():
    def __init__(self, airbox_size=300):
        self.airbox = {
                        'position_center': [0,0,0],
                        'dimensions': [airbox_size]*3
                    }
        self.coil = {
                    'position_center': [0,0,0],
                    'outer':{
                            'length': 200, # size in x direction
                            'width': 200, # size in y direction
                            'corner_radius': 50
                            },
                    'inner':{
                            'length': 150,
                            'width': 150,
                            'corner_radius': 25
                            },
                    'height': 100,
                    'mesh_size': 2,
                    }
        self.steel_channels_mesh_size = 2
        self.steel_channels = [
                                [
                                    {
                                        'position_center': [63.7, 40, 61.6],
                                        'dimensions': [123.2, 50, 3.2]
                                    },
                                    {
                                        'position_center': [63.7, 40, -61.6],
                                        'dimensions': [123.2, 50, 3.2]
                                    },
                                    {
                                        'position_center': [125.3, 40, 0],
                                        'dimensions': [3.2, 50, 126.4]
                                    }
                                ],
                                [
                                    {
                                        'position_center': [-63.7, -40, 61.6],
                                        'dimensions': [123.2, 50, 3.2]
                                    },
                                    {
                                        'position_center': [-63.7, -40, -61.6],
                                        'dimensions': [123.2, 50, 3.2]
                                    },
                                    {
                                        'position_center': [-125.3, -40, 0],
                                        'dimensions': [3.2, 50, 126.4]
                                    }
                                ],
                                [
                                    {
                                        'position_center': [0,0,0],
                                        'dimensions': [3.2, 50, 126.4]
                                    }
                                ]
                            ]
        self.search_coils = [
            *[{'position_center': [0.8, 0, z], 'dimensions': [1.6, 50, 0]} for z in [0,10,20,30,40,50,60]],
            *[{'position_center': [x, 40, 61.6], 'dimensions': [0, 50, 3.2]} for x in [2.1, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 80.0, 100.0, 110.0, 122.1]],
            *[{'position_center': [123.7, 40, z], 'dimensions': [1.6, 50, 0]} for z in [60,50,40,30,20,10,0]]
            ]

    def get_coil_rectangle_dimensions(self, type):
        x_min = self.coil['position_center'][0] - self.coil[type]['length'] / 2
        y_min = self.coil['position_center'][1] - self.coil[type]['width'] / 2
        z_min = self.coil['position_center'][1] + self.coil['height'] / 2
        return [x_min, y_min, z_min], self.coil[type]['length'], self.coil[type]['width'], self.coil[type]['corner_radius']

    def get_box_dimensions(self, box):
        x_min = box['position_center'][0] - box['dimensions'][0] / 2
        y_min = box['position_center'][1] - box['dimensions'][1] / 2
        z_min = box['position_center'][2] - box['dimensions'][2] / 2
        return [x_min, y_min, z_min], box['dimensions']

    def get_point(self, search_coil, dim_1, dim_2):
        search_coil_normal_dimension_index = search_coil['dimensions'].index(0)
        dimensions = [0,1,2]
        del dimensions[search_coil_normal_dimension_index]
        position = search_coil['position_center'].copy()
        position[dimensions[0]] += dim_1 * (search_coil['dimensions'][dimensions[0]] / 2)
        position[dimensions[1]] += dim_2 * (search_coil['dimensions'][dimensions[1]] / 2)
        return position

    def get_surface(self, geom, search_coil):
        p1 = geom.add_point(self.get_point(search_coil, +1, +1))
        p2 = geom.add_point(self.get_point(search_coil, +1, -1))
        p3 = geom.add_point(self.get_point(search_coil, -1, -1))
        p4 = geom.add_point(self.get_point(search_coil, -1, +1))
        l1 = geom.add_line(p1, p2)
        l2 = geom.add_line(p2, p3)
        l3 = geom.add_line(p3, p4)
        l4 = geom.add_line(p4, p1)
        loop = geom.add_curve_loop([l1, l2, l3, l4])
        return geom, geom.add_plane_surface(loop)

    def generate(self):
        with pygmsh.occ.Geometry() as geom:
            airbox = geom.add_box(*self.get_box_dimensions(self.airbox))
            coil_outer_boundary_surface = geom.add_rectangle(*self.get_coil_rectangle_dimensions('outer'), mesh_size=self.coil['mesh_size'])
            coil_inner_boundary_surface = geom.add_rectangle(*self.get_coil_rectangle_dimensions('inner'), mesh_size=self.coil['mesh_size'])
            flat = geom.boolean_difference(coil_outer_boundary_surface, coil_inner_boundary_surface)
            top, coil, lateral = geom.extrude(flat, [0,0,-self.coil['height']])
            geom.add_physical(coil, label="Coil")

            airbox = geom.boolean_difference(airbox, [coil], delete_other=False)
            steel_channel_volumes = []
            for i, channels in enumerate(self.steel_channels):
                channel_volume = []
                for channel in channels:
                    channel_volume.append(geom.add_box(*self.get_box_dimensions(channel)))
                if len(channel_volume) > 1:
                    steel_channel_volumes.append(geom.boolean_union(channel_volume))
                else:
                    steel_channel_volumes.append(*channel_volume)
                geom.add_physical(steel_channel_volumes[-1], label="Steel_Channel_{}".format(i))
                airbox = geom.boolean_difference(airbox, steel_channel_volumes[-1], delete_other=False)

            geom.add_physical(airbox, label="Airbox")
            search_coils = []
            for i, search_coil in enumerate(self.search_coils,1):
                geom, surface = self.get_surface(geom, search_coil)
                search_coils.append(surface)
                geom.add_physical(surface, label="Search_coil_{}".format(i))
            mesh = geom.generate_mesh()
            print(mesh)

            print("mesh.cell_data_dict", mesh.cell_data_dict)
        mesh.write("physical.vtu")

if __name__ == "__main__":
    GenerateMesh().generate()
