#!/bin/bash

# This script generate N Tetrahedron beam meshes meshes with incrising DOFs  
# The beam size is 15x15x80 and can be customized, by changing values of first and second beam corners
# The .vtu mesh files are stored in the folder ./refined_meshes/p1 or p2
# Gmsh version 4.4.1
# Meshio version 5.0.4
# Python version 3.8.10

# Execution: 
#   $./meshes_generator.bash  "cylinder/beam" 'p1'            for linear tetrahedron cylinder/beam meshes 
#   $./meshes_generator.bash  "cylinder/beam" 'q1'            for linear hexahedron cylinder/beam meshes 

if [ $1 == 'beam' ]
then 
    for i in {0..5}
    do
        # The beam first corner in world coordinates 
        declare first_beam_corner_x=-7.5
        declare first_beam_corner_y=-7.5
        declare first_beam_corner_z=0
        # The beam second corner in world coordinates 
        declare second_beam_corner_x=7.5
        declare second_beam_corner_y=7.5
        declare second_beam_corner_z=80
        # Number of pints in each direction: this is where we refine the mesh at each iteration   
        declare -i number_point_x=3+$i
        declare -i number_point_y=3+$i
        declare -i number_point_z=$number_point_x*3
        declare -i number_of_dofs=$number_point_x*$number_point_y*$number_point_z*3    
        # The file name is to be updated at each iteration to avoid erasing the previous stored mesh  
        declare file_name="./refined_meshes/beam/$2/beam_$2_$number_point_z.vtu"
        # Message prin to the user 
        echo "Creating beam $2 elements mesh  with $number_of_dofs DOFs and storing in  $file_name" 
        # Creating and sotring the mesh using the python script create_beam_mesh.py
        python3 create_beam_mesh.py --type "$2" -n $number_point_x $number_point_y $number_point_z --p0 $first_beam_corner_x $first_beam_corner_y $first_beam_corner_z --p1 $second_beam_corner_x $second_beam_corner_y $second_beam_corner_z -v -o $file_name
    done
elif [ $1 == 'cylinder' ]
then 

    declare  rad=7.5
    declare  len=80
    for size_element in  10.0 5.0 4.0 3.0 2.0 1.0 0.8 0.6
    do

        # The file name is to be updated at each iteration to avoid erasing the previous stored mesh  
        declare file_name="./refined_meshes/cylinder/$2/cylinder_$2_$size_element.vtu"
        # Message prin to the user 
        echo "Creating cylinder $2 elements mesh  with element size $size_element  and storing in  $file_name" 
        # Creating and sotring the mesh using the python script create_cylinder_mesh.py
        echo "create_cylinder_mesh.py --type "$2" --type "$2" -s $size_element -r $rad -l $len -v -o $file_name"
        python3 create_cylinder_mesh.py --type "$2" --type "$2" -s $size_element -r $rad -l $len -v -o $file_name
    done
else 
    echo "The argument $1 is invalid please choose \"cylinder\" or  \"beam\" "
fi 

