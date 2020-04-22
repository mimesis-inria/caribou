Point(1) = {-17.0648, -24.0934, 25.6331};
Point(2) = {-8.90826, -21.8884, 20.2844};
Point(3) = {-12.813,  -12.9683, 18.0072};
Point(4) = {-20.9696, -15.1734, 23.3558};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

Line Loop(5) = {2, 3, 4, 1};
Plane Surface(6) = {5};

Transfinite Surface {6} = {2, 3, 4, 1};
Transfinite Line {1, 3} = 5 Using Progression 1;
Transfinite Line {2, 4} = 5 Using Progression 1;
Reverse Surface {6};

Recombine Surface {6};

ex[] = Extrude {-4.26896, -3.94592, -8.13669} {Surface{6}; Layers{4}; Recombine;};

Physical Surface("base") = {6};
Physical Surface("top") = {ex[0]};
Physical Surface ("surface") = {6, ex[0], ex[2], ex[3], ex[4], ex[5]};
Physical Volume ("volume") = {ex[1]};

Mesh.ElementOrder = 1;
Mesh.SecondOrderLinear = 0;