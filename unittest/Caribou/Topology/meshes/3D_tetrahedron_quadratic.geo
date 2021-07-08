SetFactory("OpenCASCADE");
Sphere(1) = {5, 5, -5, 5, -Pi/2, Pi/2, 2*Pi};

Physical Surface("sphere_surface") = {1};
Physical Volume("sphere_volume") = {1};

Mesh.ElementOrder = 2;
Mesh.SecondOrderLinear = 0;