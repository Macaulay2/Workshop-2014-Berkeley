(noncommify)
(elimorder)
(setweights 1 1 3 3 2)
(setmaxdeg 6)
(simple)
vars x,y,a,b,c,h;
x*y*x-a,y*x*y-b,x*y-c;

x*y*x-a*h^2,y*x*y-b*h^2,x*y-c*h^2,
h*x-x*h,h*y-y*h,h*a-a*h,h*b-b*h,h*c-c*h;
