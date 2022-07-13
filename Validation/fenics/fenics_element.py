import basix

el = basix.create_element(basix.ElementFamily.serendipity, basix.CellType.hexahedron, 2)
res = el.points.transpose() @ el.interpolation_matrix.transpose()[:, ]
res = res.transpose()
print(res)
