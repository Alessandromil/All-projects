function Phi = Pendulum_ExpPendulum_fillphi(q,Par)

% Author  : CEIT
% Date  : 07-Nov-2022
% Time  : 16:2
% Model : Pendulum
% Version: 2.0 CEIT

Bar1_P1x = Par(1);
Bar1_P1y = Par(2);
Bar1_P1z = Par(3);
Bar1_P2x = Par(4);
Bar1_P2y = Par(5);
Bar1_P2z = Par(6);
Bar1_B1_M1x = Par(7);
Bar1_B1_M1y = Par(8);
Bar1_B1_M1z = Par(9);
Bar1_B1_M2x = Par(10);
Bar1_B1_M2y = Par(11);
Bar1_B1_M2z = Par(12);
Bar1_B1_M3x = Par(13);
Bar1_B1_M3y = Par(14);
Bar1_B1_M3z = Par(15);
Bar1_Yb1x = Par(16);
Bar1_Yb1y = Par(17);
Bar1_Yb1z = Par(18);
Bar1_ZGx = Par(19);
Bar1_ZGy = Par(20);
Bar1_ZGz = Par(21);


P2x = q(1);
P2y = q(2);
P2z = q(3);
Yb1x = q(4);
Yb1y = q(5);
Yb1z = q(6);
B1_M1x = q(7);
B1_M1y = q(8);
B1_M1z = q(9);
B1_M2x = q(10);
B1_M2y = q(11);
B1_M2z = q(12);
B1_M3x = q(13);
B1_M3y = q(14);
B1_M3z = q(15);

P1x = -0.254;
P1y = 0.912;
P1z = 0.668;
XGx = 1;
XGy = 0;
XGz = 0;
ZGx = 0;
ZGy = 0;
ZGz = 1;
YGx = 0;
YGy = 1;
YGz = 0;

Phi = zeros(15,1);
Phi(1) = Yb1x^2 + Yb1y^2 + Yb1z^2 - 1;
Phi(2) = ZGx^2 + ZGy^2 + ZGz^2 - 1;
Phi(3) = (P1x - P2x)^2 + (P1y - P2y)^2 + (P1z - P2z)^2 - 7102572929130481/576460752303423488;
Phi(4) = Yb1x*ZGx + Yb1y*ZGy + Yb1z*ZGz;
Phi(5) = - Yb1x*(P1x - P2x) - Yb1y*(P1y - P2y) - Yb1z*(P1z - P2z);
Phi(6) = - ZGx*(P1x - P2x) - ZGy*(P1y - P2y) - ZGz*(P1z - P2z);
Phi(7) = B1_M1x - (28*P1x)/37 - (9*P2x)/37 - (3*Yb1x)/500 - (9*ZGx)/500;
Phi(8) = B1_M1y - (28*P1y)/37 - (9*P2y)/37 - (3*Yb1y)/500 - (9*ZGy)/500;
Phi(9) = B1_M1z - (28*P1z)/37 - (9*P2z)/37 - (3*Yb1z)/500 - (9*ZGz)/500;
Phi(10) = B1_M2x - P2x + (17*Yb1x)/1000;
Phi(11) = B1_M2y - P2y + (17*Yb1y)/1000;
Phi(12) = B1_M2z - P2z + (17*Yb1z)/1000;
Phi(13) = B1_M3x - P2x - (13*Yb1x)/200;
Phi(14) = B1_M3y - P2y - (13*Yb1y)/200;
Phi(15) = B1_M3z - P2z - (13*Yb1z)/200;

