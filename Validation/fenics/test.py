from Caribou.Geometry import Hexahedron
import Caribou

t = Hexahedron(Caribou.Quadratic)
# for i in t.nodes():
#     print(t.L(i))
print(t.dL(t.nodes()[0]))