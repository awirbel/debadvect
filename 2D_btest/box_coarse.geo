cl1 = 5.85;
Point(1) = {0, 0, 0, cl1};
Point(2) = {0, 100, 0, cl1};
Point(3) = {100, 0, 0, cl1};
Point(4) = {100, 100, 0, cl1};
Line(1) = {2, 1};
Line(2) = {1, 3};
Line(3) = {3, 4};
Line(4) = {4, 2};
Line Loop(5) = {4, 1, 2, 3};
Plane Surface(6) = {5};
Physical Line(0) = {1, 4, 3};
Physical Line(1) = {2};
Physical Surface(9) = {6};
