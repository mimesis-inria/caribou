
def test1():
    from Caribou.Base import serialize, deserialize
    from Caribou.Mesh import cylinder

    cyl = cylinder(center1=[0, 0, 0], center2=[0, 0, 80], radius=7.5, size=24, dimension=2, quads=False)
    exported = serialize(cyl)
    assert exported == serialize(deserialize(exported))
    print "TEST PASSED"


if __name__ == '__main__':
    test1()