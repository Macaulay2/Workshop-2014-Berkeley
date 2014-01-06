(noncommify)
(setmaxdeg 8)
(setq nmodgen 3)
(algforminput)
vars x,y,z,w,c1,c2,c3;
x*y+y*x,
x*z+z*x,
y*z+z*y,
-w*y+x*w,
-w*z+y*w,
-w*x+z*w,
w^2,
c2*x+c1*y,
c3*y-c1*w,
c3*z-c2*w,
c3*w;
