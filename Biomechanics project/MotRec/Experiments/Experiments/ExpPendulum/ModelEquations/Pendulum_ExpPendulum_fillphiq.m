function Phiq = Pendulum_ExpPendulum_fillphiq(q,Par)

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
i = [ 1 1 1 3 3 3 4 4 4 5 5 5 5 5 5 6 6 6 7 7 7 8 8 8 9 9 9 10 10 10 11 11 11 12 12 12 13 13 13 14 14 14 15 15 15 ];

j = [ 4 5 6 1 2 3 4 5 6 1 2 3 4 5 6 1 2 3 1 4 7 2 5 8 3 6 9 1 4 10 2 5 11 3 6 12 1 4 13 2 5 14 3 6 15 ];

s = [ 2*Yb1x
2*Yb1y
2*Yb1z
2*P2x - 2*P1x
2*P2y - 2*P1y
2*P2z - 2*P1z
ZGx
ZGy
ZGz
Yb1x
Yb1y
Yb1z
P2x - P1x
P2y - P1y
P2z - P1z
ZGx
ZGy
ZGz
-9/37
-3/500
1
-9/37
-3/500
1
-9/37
-3/500
1
-1
17/1000
1
-1
17/1000
1
-1
17/1000
1
-1
-13/200
1
-1
-13/200
1
-1
-13/200
1
];

Phiq = sparse(i,j,s,15,15);

