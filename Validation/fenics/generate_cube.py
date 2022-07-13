# Import modules:
import gmsh

element = "hexa"
order = 2
GUI = False

# Initialize gmsh:
gmsh.initialize()

# cube points:
L, W, H = 80, 15, 15
lc = L/6

if element == "hexa":
    point1 = gmsh.model.geo.add_point(-H / 2, -W / 2, 0.0, lc)
    point2 = gmsh.model.geo.add_point(H / 2, -W / 2, 0.0, lc)
    point3 = gmsh.model.geo.add_point(H / 2, W / 2, 0.0, lc)
    point4 = gmsh.model.geo.add_point(-H / 2, W / 2, 0, lc)

    # Edge of cube:
    line1 = gmsh.model.geo.add_line(point1, point2)
    line2 = gmsh.model.geo.add_line(point2, point3)
    line3 = gmsh.model.geo.add_line(point3, point4)
    line4 = gmsh.model.geo.add_line(point4, point1)

    # faces of cube:
    face1 = gmsh.model.geo.add_curve_loop([line1, line2, line3, line4])

    # surfaces of cube:
    s1 = gmsh.model.geo.add_plane_surface([face1])
    gmsh.model.geo.synchronize()

    gmsh.option.setNumber("Mesh.RecombineAll", 1)
    gmsh.option.setNumber("Mesh.Algorithm", 8)

    sl = gmsh.model.geo.addSurfaceLoop([s1])

    gmsh.model.geo.synchronize()

    # Create the relevant Gmsh data structures
    # from Gmsh model.
    gmsh.option.setNumber("Mesh.ElementOrder", order)
    if order == 2:
        gmsh.option.setNumber("Mesh.SecondOrderIncomplete", 1)

    # Generate mesh:
    gmsh.model.mesh.generate(2)
    ov = gmsh.model.geo.extrude([(2, 1)], 0, 0, L, [lc/2, lc/2], [0.5, 1], recombine=True)
    gmsh.model.geo.synchronize()
    gmsh.model.mesh.generate(3)

    gmsh.model.addPhysicalGroup(3, [1], 187)
    gmsh.model.setPhysicalName(3, 187, 'Volume')
    gmsh.model.addPhysicalGroup(2, [1, 13, 17, 21, 25, 26], 188)
    gmsh.model.setPhysicalName(2, 188, 'Surface')


elif element == "tetra":
    # cube points:
    point1 = gmsh.model.geo.add_point(-H / 2, -W / 2, 0.0, lc)
    point2 = gmsh.model.geo.add_point(H / 2, -W / 2, 0.0, lc)
    point3 = gmsh.model.geo.add_point(H / 2, W / 2, 0.0, lc)
    point4 = gmsh.model.geo.add_point(-H / 2, W / 2, 0, lc)
    point5 = gmsh.model.geo.add_point(-H / 2, W / 2, L, lc)
    point6 = gmsh.model.geo.add_point(-H / 2, -W / 2, L, lc)
    point7 = gmsh.model.geo.add_point(H / 2, -W / 2, L, lc)
    point8 = gmsh.model.geo.add_point(H / 2, W / 2, L, lc)

    # Edge of cube:
    line1 = gmsh.model.geo.add_line(point1, point2)
    line2 = gmsh.model.geo.add_line(point2, point3)
    line3 = gmsh.model.geo.add_line(point3, point4)
    line4 = gmsh.model.geo.add_line(point4, point1)
    line5 = gmsh.model.geo.add_line(point5, point6)
    line6 = gmsh.model.geo.add_line(point6, point7)
    line7 = gmsh.model.geo.add_line(point7, point8)
    line8 = gmsh.model.geo.add_line(point8, point5)
    line9 = gmsh.model.geo.add_line(point4, point5)
    line10 = gmsh.model.geo.add_line(point6, point1)
    line11 = gmsh.model.geo.add_line(point7, point2)
    line12 = gmsh.model.geo.add_line(point3, point8)

    # faces of cube:
    face1 = gmsh.model.geo.add_curve_loop([line1, line2, line3, line4])
    face2 = gmsh.model.geo.add_curve_loop([line5, line6, line7, line8])
    face3 = gmsh.model.geo.add_curve_loop([line9, line5, line10, -line4])
    face4 = gmsh.model.geo.add_curve_loop([line9, -line8, -line12, line3])
    face5 = gmsh.model.geo.add_curve_loop([line6, line11, -line1, -line10])
    face6 = gmsh.model.geo.add_curve_loop([line11, line2, line12, -line7])

    # surfaces of cube:
    s1 = gmsh.model.geo.add_plane_surface([face1])
    s2 = gmsh.model.geo.add_plane_surface([face2])
    s3 = gmsh.model.geo.add_plane_surface([face3])
    s4 = gmsh.model.geo.add_plane_surface([face4])
    s5 = gmsh.model.geo.add_plane_surface([face5])
    s6 = gmsh.model.geo.add_plane_surface([face6])
    gmsh.model.geo.synchronize()

    shells = []

    sl = gmsh.model.geo.addSurfaceLoop([s1, s2, s3, s4, s5, s6])
    shells.append(sl)

    gmsh.model.geo.addVolume(shells, 186)

    gmsh.model.addPhysicalGroup(3, [186], 187)
    gmsh.model.setPhysicalName(3, 187, 'Volume')
    gmsh.model.addPhysicalGroup(2, [s1, s2, s3, s4, s5, s6], 188)
    gmsh.model.setPhysicalName(2, 188, 'Surface')
    gmsh.model.addPhysicalGroup(2, [s1], 1)
    gmsh.model.setPhysicalName(2, 188, 'left')

    gmsh.model.geo.synchronize()

    # Create the relevant Gmsh data structures
    # from Gmsh model.
    gmsh.option.setNumber("Mesh.ElementOrder", order)
    # gmsh.option.

    # Generate mesh:
    gmsh.model.mesh.generate(3)

else:
    assert "Unknown element"

# Write mesh data:
gmsh.write("../meshes/fenics/cube" + "_" + element + "_" + str(order) + ".msh")

# Creates  graphical user interface
if GUI:
    gmsh.fltk.run()

# It finalize the Gmsh API
gmsh.finalize()
