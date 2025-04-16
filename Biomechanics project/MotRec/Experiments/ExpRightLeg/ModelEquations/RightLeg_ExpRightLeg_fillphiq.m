function Phiq = RightLeg_ExpRightLeg_fillphiq(q,Par)

% Author  : CEIT
% Date  : 29-Nov-2024
% Time  : 11:44
% Model : RightLeg
% Version: 2.0 CEIT

Shank_AJCx = Par(1);
Shank_AJCy = Par(2);
Shank_AJCz = Par(3);
Shank_KJCx = Par(4);
Shank_KJCy = Par(5);
Shank_KJCz = Par(6);
Shank_MarkerSet_G10_RSx = Par(7);
Shank_MarkerSet_G10_RSy = Par(8);
Shank_MarkerSet_G10_RSz = Par(9);
Shank_MarkerSet_G10_LMx = Par(10);
Shank_MarkerSet_G10_LMy = Par(11);
Shank_MarkerSet_G10_LMz = Par(12);
Shank_MarkerSet_G10_MMx = Par(13);
Shank_MarkerSet_G10_MMy = Par(14);
Shank_MarkerSet_G10_MMz = Par(15);
Shank_Xsx = Par(16);
Shank_Xsy = Par(17);
Shank_Xsz = Par(18);
Shank_Zsx = Par(19);
Shank_Zsy = Par(20);
Shank_Zsz = Par(21);
Foot_AJCx = Par(22);
Foot_AJCy = Par(23);
Foot_AJCz = Par(24);
Foot_MarkerSet_G10_F2x = Par(25);
Foot_MarkerSet_G10_F2y = Par(26);
Foot_MarkerSet_G10_F2z = Par(27);
Foot_MarkerSet_G10_F1x = Par(28);
Foot_MarkerSet_G10_F1y = Par(29);
Foot_MarkerSet_G10_F1z = Par(30);
Foot_MarkerSet_G10_F3x = Par(31);
Foot_MarkerSet_G10_F3y = Par(32);
Foot_MarkerSet_G10_F3z = Par(33);
Foot_Xfx = Par(34);
Foot_Xfy = Par(35);
Foot_Xfz = Par(36);
Foot_Yfx = Par(37);
Foot_Yfy = Par(38);
Foot_Yfz = Par(39);
Foot_Zfx = Par(40);
Foot_Zfy = Par(41);
Foot_Zfz = Par(42);
Thigh_KJCx = Par(43);
Thigh_KJCy = Par(44);
Thigh_KJCz = Par(45);
Thigh_HJCx = Par(46);
Thigh_HJCy = Par(47);
Thigh_HJCz = Par(48);
Thigh_MarkerSet_G10_RTx = Par(49);
Thigh_MarkerSet_G10_RTy = Par(50);
Thigh_MarkerSet_G10_RTz = Par(51);
Thigh_MarkerSet_G10_LFEx = Par(52);
Thigh_MarkerSet_G10_LFEy = Par(53);
Thigh_MarkerSet_G10_LFEz = Par(54);
Thigh_MarkerSet_G10_MFEx = Par(55);
Thigh_MarkerSet_G10_MFEy = Par(56);
Thigh_MarkerSet_G10_MFEz = Par(57);
Thigh_Xtx = Par(58);
Thigh_Xty = Par(59);
Thigh_Xtz = Par(60);
Thigh_Ztx = Par(61);
Thigh_Zty = Par(62);
Thigh_Ztz = Par(63);

AJCx = q(1);
AJCy = q(2);
AJCz = q(3);
KJCx = q(4);
KJCy = q(5);
KJCz = q(6);
HJCx = q(7);
HJCy = q(8);
HJCz = q(9);
Xsx = q(10);
Xsy = q(11);
Xsz = q(12);
Zsx = q(13);
Zsy = q(14);
Zsz = q(15);
Xfx = q(16);
Xfy = q(17);
Xfz = q(18);
Yfx = q(19);
Yfy = q(20);
Yfz = q(21);
Zfx = q(22);
Zfy = q(23);
Zfz = q(24);
Xtx = q(25);
Xty = q(26);
Xtz = q(27);
Ztx = q(28);
Zty = q(29);
Ztz = q(30);
MarkerSet_G10_RSx = q(31);
MarkerSet_G10_RSy = q(32);
MarkerSet_G10_RSz = q(33);
MarkerSet_G10_LMx = q(34);
MarkerSet_G10_LMy = q(35);
MarkerSet_G10_LMz = q(36);
MarkerSet_G10_MMx = q(37);
MarkerSet_G10_MMy = q(38);
MarkerSet_G10_MMz = q(39);
MarkerSet_G10_F2x = q(40);
MarkerSet_G10_F2y = q(41);
MarkerSet_G10_F2z = q(42);
MarkerSet_G10_F1x = q(43);
MarkerSet_G10_F1y = q(44);
MarkerSet_G10_F1z = q(45);
MarkerSet_G10_F3x = q(46);
MarkerSet_G10_F3y = q(47);
MarkerSet_G10_F3z = q(48);
MarkerSet_G10_RTx = q(49);
MarkerSet_G10_RTy = q(50);
MarkerSet_G10_RTz = q(51);
MarkerSet_G10_LFEx = q(52);
MarkerSet_G10_LFEy = q(53);
MarkerSet_G10_LFEz = q(54);
MarkerSet_G10_MFEx = q(55);
MarkerSet_G10_MFEy = q(56);
MarkerSet_G10_MFEz = q(57);
OGx = 0;
OGy = 0;
OGz = 0;
XGx = 1;
XGy = 0;
XGz = 0;
ZGx = 0;
ZGy = 0;
ZGz = 1;
YGx = 0;
YGy = 1;
YGz = 0;
i = [ 1 1 1 2 2 2 3 3 3 3 3 3 4 4 4 4 4 4 5 5 5 5 5 5 5 5 5 6 6 6 6 6 6 6 6 6 7 7 7 7 7 8 8 8 8 8 9 9 9 9 9 10 10 10 10 10 11 11 11 11 11 12 12 12 12 12 13 13 13 13 13 14 14 14 14 14 15 15 15 15 15 16 16 16 17 17 17 18 18 18 19 19 19 19 19 19 20 20 20 20 20 20 21 21 21 21 21 21 22 22 22 22 22 23 23 23 23 23 24 24 24 24 24 25 25 25 25 25 26 26 26 26 26 27 27 27 27 27 28 28 28 28 28 29 29 29 29 29 30 30 30 30 30 31 31 31 32 32 32 33 33 33 33 33 33 34 34 34 34 34 34 35 35 35 35 35 35 35 35 35 36 36 36 36 36 36 36 36 36 37 37 37 37 37 38 38 38 38 38 39 39 39 39 39 40 40 40 40 41 41 41 41 42 42 42 42 43 43 43 43 44 44 44 44 45 45 45 45 46 46 46 46 46 46 47 47 48 48 49 49 ];

j = [ 10 11 12 13 14 15 1 2 3 4 5 6 10 11 12 13 14 15 1 2 3 4 5 6 10 11 12 1 2 3 4 5 6 13 14 15 1 4 10 13 31 2 5 11 14 32 3 6 12 15 33 1 4 10 13 34 2 5 11 14 35 3 6 12 15 36 1 4 10 13 37 2 5 11 14 38 3 6 12 15 39 16 17 18 19 20 21 22 23 24 16 17 18 19 20 21 16 17 18 22 23 24 19 20 21 22 23 24 1 16 19 22 40 2 17 20 23 41 3 18 21 24 42 1 16 19 22 43 2 17 20 23 44 3 18 21 24 45 1 16 19 22 46 2 17 20 23 47 3 18 21 24 48 25 26 27 28 29 30 4 5 6 7 8 9 25 26 27 28 29 30 4 5 6 7 8 9 25 26 27 4 5 6 7 8 9 28 29 30 4 7 25 28 49 5 8 26 29 50 6 9 27 30 51 4 7 28 52 5 8 29 53 6 9 30 54 4 7 28 55 5 8 29 56 6 9 30 57 13 14 15 16 17 18 13 28 14 29 15 30 ];

s = [ 2*Xsx
2*Xsy
2*Xsz
2*Zsx
2*Zsy
2*Zsz
2*AJCx - 2*KJCx
2*AJCy - 2*KJCy
2*AJCz - 2*KJCz
2*KJCx - 2*AJCx
2*KJCy - 2*AJCy
2*KJCz - 2*AJCz
Zsx
Zsy
Zsz
Xsx
Xsy
Xsz
-Xsx
-Xsy
-Xsz
Xsx
Xsy
Xsz
KJCx - AJCx
KJCy - AJCy
KJCz - AJCz
-Zsx
-Zsy
-Zsz
Zsx
Zsy
Zsz
KJCx - AJCx
KJCy - AJCy
KJCz - AJCz
-2445/3658
-1213/3658
-301/10000
-157/10000
1
-2445/3658
-1213/3658
-301/10000
-157/10000
1
-2445/3658
-1213/3658
-301/10000
-157/10000
1
-3749/3658
91/3658
33/2000
-193/5000
1
-3749/3658
91/3658
33/2000
-193/5000
1
-3749/3658
91/3658
33/2000
-193/5000
1
-1784/1829
-45/1829
-33/2000
193/5000
1
-1784/1829
-45/1829
-33/2000
193/5000
1
-1784/1829
-45/1829
-33/2000
193/5000
1
2*Xfx
2*Xfy
2*Xfz
2*Yfx
2*Yfy
2*Yfz
2*Zfx
2*Zfy
2*Zfz
Yfx
Yfy
Yfz
Xfx
Xfy
Xfz
Zfx
Zfy
Zfz
Xfx
Xfy
Xfz
Zfx
Zfy
Zfz
Yfx
Yfy
Yfz
-1
-611/10000
-29/10000
-13/400
1
-1
-611/10000
-29/10000
-13/400
1
-1
-611/10000
-29/10000
-13/400
1
-1
-1307/10000
47/1250
59/10000
1
-1
-1307/10000
47/1250
59/10000
1
-1
-1307/10000
47/1250
59/10000
1
-1
-833/10000
41/1000
-429/5000
1
-1
-833/10000
41/1000
-429/5000
1
-1
-833/10000
41/1000
-429/5000
1
2*Xtx
2*Xty
2*Xtz
2*Ztx
2*Zty
2*Ztz
2*KJCx - 2*HJCx
2*KJCy - 2*HJCy
2*KJCz - 2*HJCz
2*HJCx - 2*KJCx
2*HJCy - 2*KJCy
2*HJCz - 2*KJCz
Ztx
Zty
Ztz
Xtx
Xty
Xtz
-Xtx
-Xty
-Xtz
Xtx
Xty
Xtz
HJCx - KJCx
HJCy - KJCy
HJCz - KJCz
-Ztx
-Zty
-Ztz
Ztx
Zty
Ztz
HJCx - KJCx
HJCy - KJCy
HJCz - KJCz
-3453/4700
-1247/4700
-491/5000
3/500
1
-3453/4700
-1247/4700
-491/5000
3/500
1
-3453/4700
-1247/4700
-491/5000
3/500
1
-232/235
-3/235
-11/200
1
-232/235
-3/235
-11/200
1
-232/235
-3/235
-11/200
1
-238/235
3/235
11/200
1
-238/235
3/235
11/200
1
-238/235
3/235
11/200
1
Xfx
Xfy
Xfz
Zsx
Zsy
Zsz
-1
1
-1
1
-1
1
];

Phiq = sparse(i,j,s,49,57);

