function rQ = DSolver(a,b,c)
b3a = b./(3*a);
ca = c./a;
b3a2 = b3a.^2;
R0 = b3a2.*b3a - b3a.*ca./2;
Q = b3a2 - ca./3;
rQ = sqrt(Q);
rQ3 = sqrt(Q.^3);

% R = R0 + .5*d/a;
% th = arccos(R/rQ3) + 2*pi;
% x2 = -2*rQ*cos(th/3) - b3a;
end