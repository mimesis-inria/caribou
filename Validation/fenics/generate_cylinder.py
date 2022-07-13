# Import modules:
import gmsh

element = "tetra"
order = 1
GUI = True

# Initialize gmsh:
gmsh.initialize()

# cube points:
L, W, H = 80, 15, 15
lc = L / 20

if element == "hexa":
    point1 = gmsh.model.geo.add_point(-H / 2, 0, 0.0, lc)
    point2 = gmsh.model.geo.add_point(0, -W / 2, 0.0, lc)
    point3 = gmsh.model.geo.add_point(H / 2, 0, 0.0, lc)
    point4 = gmsh.model.geo.add_point(0, W / 2, 0, lc)
    point5 = gmsh.model.geo.add_point(0, 0, 0, lc)

    # Edge of cube:
    line1 = gmsh.model.geo.add_circle_arc(point1, point5, point2)
    line2 = gmsh.model.geo.add_circle_arc(point2, point5, point3)
    line3 = gmsh.model.geo.add_circle_arc(point3, point5, point4)
    line4 = gmsh.model.geo.add_circle_arc(point4, point5, point1)

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
    ov = gmsh.model.geo.extrude([(2, 1)], 0, 0, L, [lc / 2, lc / 2], [0.5, 1], recombine=True)
    gmsh.model.geo.synchronize()
    gmsh.model.mesh.generate(3)

    gmsh.model.addPhysicalGroup(3, [1], 187)
    gmsh.model.setPhysicalName(3, 187, 'Volume')
    gmsh.model.addPhysicalGroup(2, [1, 13, 17, 21, 25, 26], 188)
    gmsh.model.setPhysicalName(2, 188, 'Surface')


elif element == "tetra":
    gmsh.model.occ.addCylinder(0, 0, 0, 0, 0, L, W / 2)
    gmsh.model.occ.synchronize()

    # Create the relevant Gmsh data structures
    # from Gmsh model.
    gmsh.option.setNumber("Mesh.ElementOrder", order)
    # gmsh.option.

    # Generate mesh:
    gmsh.model.mesh.generate(3)
    gmsh.model.addPhysicalGroup(3, [1], 187)
    gmsh.model.setPhysicalName(3, 187, 'Volume')
    gmsh.model.addPhysicalGroup(2, [1, 2, 3], 188)
    gmsh.model.setPhysicalName(2, 188, 'Surface')

else:
    assert "Unknown element"

# Write mesh data:
gmsh.write("../meshes/fenics/cylinder" + "_" + element + "_" + str(order) + ".msh")

# Creates  graphical user interface
if GUI:
    gmsh.fltk.run()

# It finalize the Gmsh API
gmsh.finalize()
