Point(1) =  {-5, 0, 0, 1.0};
Point(2) =  {-4, 0, 0, 1.0};
Point(3) =  {-3, 0, 0, 1.0};
Point(4) =  {-2, 0, 0, 1.0};
Point(5) =  {-1, 0, 0, 1.0};
Point(6) =  { 0, 0, 0, 1.0};
Point(7) =  { 1, 0, 0, 1.0};
Point(8) =  { 2, 0, 0, 1.0};
Point(9) =  { 3, 0, 0, 1.0};
Point(10) = { 4, 0, 0, 1.0};
Point(11) = { 5, 0, 0, 1.0};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 5};
Line(5) = {5, 6};
Line(6) = {6, 7};
Line(7) = {7, 8};
Line(8) = {8, 9};
Line(9) = {9, 10};
Line(10) = {10, 11};

Physical Point("starting_point") = {1};
Physical Point("ending_point") = {11};
Physical Curve("Line") = {1,2,3,4,5,6,7,8,9,10};

Mesh.ElementOrder = 1;
Mesh.SecondOrderLinear = 1;
Mesh 1;
