function [X, Y, Z] = cylinder2P(R, N,r1,r2)

% The parametric surface will consist of a series of N-sided
% polygons with successive radii given by the array R.
% Z increases in equal sized steps from 0 to 1.
% Set up an array of angles for the polygon.

theta = linspace(0,2*pi,N);
m = length(R); % Number of radius values
% supplied.

if m == 1 % Only one radius value supplied.
    R = [R; R]; % Add a duplicate radius to make
    m = 2; % a cylinder.
end

X = zeros(m, N); % Preallocate memory.
Y = zeros(m, N);
Z = zeros(m, N);

%Transpose vectors
r1=r1';
r2=r2';

v=(r2-r1)/sqrt((r2-r1)*(r2-r1)'); %Normalized vector;
%cylinder axis described by: r(t)=r1+v*t for 0<t<1
R2=rand(1,3); %linear independent vector (of v)
x2=v-R2/(R2*v'); %orthogonal vector to v
x2=x2/sqrt(x2*x2'); %orthonormal vector to v
x3=cross(v,x2); %vector orthonormal to v and x2
x3=x3/sqrt(x3*x3');
r1x=r1(1);r1y=r1(2);r1z=r1(3);
r2x=r2(1);r2y=r2(2);r2z=r2(3);
vx=v(1);vy=v(2);vz=v(3);
x2x=x2(1);x2y=x2(2);x2z=x2(3);
x3x=x3(1);x3y=x3(2);x3z=x3(3);
time=linspace(0,1,m);
for j = 1 : m
    t=time(j);
    X(j, :) = r1x+(r2x-r1x)*t+R(j)*cos(theta)*x2x+R(j)*sin(theta)*x3x;
    Y(j, :) = r1y+(r2y-r1y)*t+R(j)*cos(theta)*x2y+R(j)*sin(theta)*x3y;
    Z(j, :) = r1z+(r2z-r1z)*t+R(j)*cos(theta)*x2z+R(j)*sin(theta)*x3z;
end
end