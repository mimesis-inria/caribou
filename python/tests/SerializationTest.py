
def test1():
    from Caribou.Base import deserialize
    from Caribou.Mesh import cylinder

    cyl = cylinder(center1=[0, 0, 0], center2=[0, 0, 80], radius=7.5, size=24, dimension=2, quads=False)
    cyl2 = deserialize(cyl.serialize())
    assert cyl.serialize() == cyl2.serialize()


if __name__ == '__main__':
    test1()