! ----------------------------------------------------
! Playback File generated automaticaly - CEIT 2003.
! ----------------------------------------------------
Simulation Playback
ObjectMotion=myMotion File=ExpPendulum_Pendulum_Pend_simple01.sim

! ----------------------------------------------------
! Dimensions of the graphic objects
! ----------------------------------------------------
! Sphere radius of the bony palpable landmarks
Parameter = LmarkRadius   Value = 0.007;
! Sphere radius of each joint centre
Parameter = JointRadius   Value = 0.015;
! Sphere radius of the skin-markers
Parameter = MarkerRadius   Value = 0.015;
!Radius of the cylinder that connects each experimental skin-marker and body skin-marker
Parameter = AuxWireRadius   Value = 0.004;
!Radius of the cylinder that represents each muscles of the model
Parameter = MuscleRadius   Value = 0.006;

! ----------------------------------------------------
!Graphics files (bodies, joints, landmarks)
! ----------------------------------------------------

! ----------------------------------------------------
!Graphic file - Body graphics 
!----------------------------------------------------

! Wireframe graphics between P1-P2 in Body Bar1
Graphic=Bar1_WireP1P2 CYLINDER x0=0  y0=0  z0=0  x1=0.111  y1=0  z1=0  Radius=0.01;  Facets=40  Material=Gray-25Object=Bar1   Entity=Bar1_WireP1P2   Motion=myMotion

! Wireframe graphics between P2-B1_M3 in Body Bar1
Graphic=Bar1_WireP2B1_M3 CYLINDER x0=0.111  y0=0  z0=0  x1=0.111  y1=0.065  z1=0  Radius=0.01;  Facets=40  Material=Gray-25Object=Bar1   Entity=Bar1_WireP2B1_M3   Motion=myMotion

! ----------------------------------------------------
!Graphic file - Joints 
!----------------------------------------------------

!---------------------------------------
! Body: Bar1
!---------------------------------------

! P1 graphics in Body Bar1
Graphic=Bar1_P1 SPHERE x0=0  y0=0  z0=0  Radius=JointRadius;  Facets=40  Material=red
Object=Bar1   Entity=Bar1_P1   Motion=myMotion

! ----------------------------------------------------
! Markers in the Body(Rigidly fixed to Bodies) and Measured 
!----------------------------------------------------

!-------------------------------------------------------------
! Graphics File:  Markers of the Bodies
!-------------------------------------------------------------

! Generic sphere for measured markers
Graphic=MeasuredMarker SPHERE x0=0.0  y0=0.0  z0=0.0  Radius=MarkerRadius;  Facets=40  Material=green

! B1_M1 graphics in Body Bar1
Graphic=Bar1_B1_M1 SPHERE x0=0.027  y0=0.006  z0=0.018  Radius=MarkerRadius;  Facets=40  Material=blue
Object=Bar1   Entity=Bar1_B1_M1   Motion=myMotion

! B1_M2 graphics in Body Bar1
Graphic=Bar1_B1_M2 SPHERE x0=0.111  y0=-0.017  z0=0  Radius=MarkerRadius;  Facets=40  Material=blue
Object=Bar1   Entity=Bar1_B1_M2   Motion=myMotion

! B1_M3 graphics in Body Bar1
Graphic=Bar1_B1_M3 SPHERE x0=0.111  y0=0.065  z0=0  Radius=MarkerRadius;  Facets=40  Material=blue
Object=Bar1   Entity=Bar1_B1_M3   Motion=myMotion


!-------------------------------------------------------------
! Graphics File:  Markers measured with VICON
!-------------------------------------------------------------

! Spheres
Object=MeasuredMarker_B1_M1  Entity=MeasuredMarker   Motion=myMotion
Object=MeasuredMarker_B1_M2  Entity=MeasuredMarker   Motion=myMotion
Object=MeasuredMarker_B1_M3  Entity=MeasuredMarker   Motion=myMotion
!Cylinder between each experimental skin marker and each body skin marker 
Graphic=Aux_Wireframe CYLINDER x0=0  y0=0  z0=0  x1=1  y1=0  z1=0  Radius=AuxWireRadius;  Facets=40  Material=red
Object=Aux_B1_M1   Entity=Aux_Wireframe   Motion=myMotion

Object=Aux_B1_M2   Entity=Aux_Wireframe   Motion=myMotion

Object=Aux_B1_M3   Entity=Aux_Wireframe   Motion=myMotion

