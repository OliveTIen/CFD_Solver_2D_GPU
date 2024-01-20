SetFactory("OpenCASCADE");

mesh_size= 0.5;
w=1.0;
r = 0.1;

Point(1) = {-w,-w,0,mesh_size};
Point(2) = {2*w,-w,0,mesh_size}; 
Point(3) = {2*w,w,0,mesh_size};
Point(4) = {-w,w,0,mesh_size};
Line(1) = {1,2} ;
Line(2) = {2,3} ;
Line(3) = {3,4} ;
Line(4) = {4,1} ;
Circle(5) = {0, 0, 0, r, 0, 2*Pi};
Curve Loop(6) = {1, 2, 3, 4};
Curve Loop(7) = {5};
Plane Surface(1) = {6, 7};

Physical Curve("inf", 8) = {1,2,3,4};
Physical Curve("obstacle", 11) = {5};
Physical Surface("inner", 12) = {1};

Transfinite Curve {1} = 61 Using Progression 1;
Transfinite Curve {2} = 41 Using Progression 1;
Transfinite Curve {3} = 61 Using Progression 1;
Transfinite Curve {4} = 41 Using Progression 1;
Transfinite Curve {5} = 41 Using Progression 1;
