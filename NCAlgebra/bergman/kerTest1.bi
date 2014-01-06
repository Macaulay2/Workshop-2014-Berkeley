(noncommify)
(setweights 1 1 1 1 2 2 2 1)
(setmaxdeg 8)
(setq nmodgen 4)
(algforminput)
vars x,y,z,w,c1,c2,c3,r1;
x*y+y*x,
x*z+z*x,
y*z+z*y,
-w*y+x*w,
-w*z+y*w,
-w*x+z*w,
w^2,
r1*x-c1,
r1*y-c2,
r1*w-c3;
