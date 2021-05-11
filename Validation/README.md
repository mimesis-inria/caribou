# Validation

## Creation of the meshes
Validation are performed on rectangular beam meshes. To generate these meshes, the following 
dependencies must be installed:
```shell
$ apt install gmsh
$ pip3 install --user numpy meshio
```

The validation meshes can thereafter be created using:
```shell
$ cd meshes

$ # Rectangular beams

$ # Linear tetrahedral mesh (4 nodes)
$ python3 create_beam_mesh.py -t p1 -n 3 3 9 --p0 -7.5 -7.5 0 --p1 7.5 7.5 80 -o beam_p1.vtu
$ # Quadratic tetrahedral mesh (10 nodes)
$ python3 create_beam_mesh.py -t p2 -n 3 3 9 --p0 -7.5 -7.5 0 --p1 7.5 7.5 80 -o beam_p2.vtu
$ # Trilinear hexahedral mesh (8 nodes)
$ python3 create_beam_mesh.py -t q1 -n 3 3 9 --p0 -7.5 -7.5 0 --p1 7.5 7.5 80 -o beam_q1.vtu
$ # Quadratic hexahedral mesh (20 nodes)
$ python3 create_beam_mesh.py -t q2 -n 3 3 9 --p0 -7.5 -7.5 0 --p1 7.5 7.5 80 -o beam_q2.vtu

$ # Cylindrical beam

$ # Linear tetrahedral mesh (4 nodes)
$ python create_cylinder_mesh.py -t p1 -s 0.18 -r 1 -l 3 -v -o cylinder_p1.vtu 
```

## Running the fenics validation scripts
Scripts starting with "fenics_[...].py" are executed using the [FEniCS Project](https://fenicsproject.org). In order to 
run these validation scenes, `podman` or `docker` must be installed. At the top of each FEniCS scripts in this directory,
the docker/podman line used to compute the solution is written.

## Running the FEBio validation scripts
Scripts starting with "febio_[...].py" are executed using [FEBio](https://febio.org/) v3.3.0. To reproduce
the validation data, the docker image can be created using the following command:
```shell
$ docker build -t febio_validation - < febio_validation.dockerfile
```

The docker/podman line used to compute the solution from this image is written at the top of each
FEBio scripts.

## Running the manufactured solution script
The manufactured solution is computed using sympy and meshio. To reproduce
the validation data, the docker image can be created using the following command:
```shell
$ docker build -t caribou_validation - < manufactured_solution/caribou.dockerfile
```

The docker/podman line used to compute the solution from this image is written at the top of the
script.
