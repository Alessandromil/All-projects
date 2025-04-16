function [X,Y,Z] = drawsphere(r,x0,y0,z0)

%Create a unit sphere centered in the origin
[x1,y1,z1] = sphere;

%Scale my multiplying by the radius
x2 = x1*r;
y2 = y1*r;
z2 = z1*r;

X = x2+x0;
Y = y2+y0;
Z = z2+z0;

end