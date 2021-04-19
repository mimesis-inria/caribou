Point(1) = {-5, -5, 2};
Point(2) = { 5, -5, 2};
Point(3) = { 5,  5, 2};
Point(4) = {-5,  5, 2};

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

Mesh.ElementOrder = 1;
Mesh.SecondOrderLinear = 0;
Mesh 2;