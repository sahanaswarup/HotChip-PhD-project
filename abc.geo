Mesh.Algorithm3D=1;
lc = 0.2;
a = 2;
b = 1;
h = 2;

Point(1) = {-0.5*a,   0, 0, lc};
Point(2) = { 0.5*a,   0, 0, lc};
Point(3) = { 0.5*a+b, 0, 0, lc};
Point(4) = { 0.5*a+b, h, 0, lc};
Point(5) = { 0.5*a,   h, 0, lc};
Point(6) = {-0.5*a,   h, 0, lc};
Point(7) = {-0.5*a-b, h, 0, lc};
Point(8) = {-0.5*a-b, 0, 0, lc};
Point(9) = {-0.5*a,   0, 0, lc};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 5};
Line(5) = {5, 6};
Line(6) = {6, 7};
Line(7) = {7, 8};
Line(8) = {8, 1};
Line(9) = {2, 5};
Line(10) = {1, 6};

Line Loop(11) = {1, 9, 5, -10};
Plane Surface(12) = {11};
Line Loop(13) = {9, -4, -3, -2};
Plane Surface(14) = {13};
Line Loop(15) = {8, 10, 6, 7};
Plane Surface(16) = {15};

out[] = Extrude {0,0,5} {
    Surface {12};
};
Color {255,255,0} {Surface{out[0], out[2], out[3], out[4], out[5]};}
Color {255,255,0} {Volume{out[1]};}
Color {255,255,0} {Surface {12};}
Color {0,255,0} {Surface {14, 16};}


