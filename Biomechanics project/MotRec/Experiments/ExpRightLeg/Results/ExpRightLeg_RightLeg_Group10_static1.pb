! ----------------------------------------------------
! Playback File generated automaticaly - CEIT 2003.
! ----------------------------------------------------
Simulation Playback
ObjectMotion=myMotion File=ExpRightLeg_RightLeg_Group10_static1.sim

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

! Wireframe graphics between AJC-KJC in Body Shank
Graphic=Shank_WireAJCKJC CYLINDER x0=0  y0=0  z0=0  x1=0  y1=0.3658  z1=0  Radius=0.01;  Facets=40  Material=Gray-25Object=Shank   Entity=Shank_WireAJCKJC   Motion=myMotion

! Wireframe graphics between AJC-MarkerSet_G10_F3 in Body Foot
Graphic=Foot_WireAJCMarkerSet_G10_F3 CYLINDER x0=0  y0=0  z0=0  x1=0.0833  y1=-0.041  z1=0.0858  Radius=0.01;  Facets=40  Material=Gray-25Object=Foot   Entity=Foot_WireAJCMarkerSet_G10_F3   Motion=myMotion

! Wireframe graphics between KJC-HJC in Body Thigh
Graphic=Thigh_WireKJCHJC CYLINDER x0=0  y0=0  z0=0  x1=0  y1=0.47  z1=0  Radius=0.01;  Facets=40  Material=Gray-25Object=Thigh   Entity=Thigh_WireKJCHJC   Motion=myMotion

! ----------------------------------------------------
!Graphic file - Joints 
!----------------------------------------------------

!---------------------------------------
! Body: Foot
!---------------------------------------

! AJC graphics in Body Foot
Graphic=Foot_AJC SPHERE x0=0  y0=0  z0=0  Radius=JointRadius;  Facets=40  Material=red
Object=Foot   Entity=Foot_AJC   Motion=myMotion

!---------------------------------------
! Body: Shank
!---------------------------------------

! KJC graphics in Body Shank
Graphic=Shank_KJC SPHERE x0=0  y0=0.3658  z0=0  Radius=JointRadius;  Facets=40  Material=red
Object=Shank   Entity=Shank_KJC   Motion=myMotion

! ----------------------------------------------------
! Markers in the Body(Rigidly fixed to Bodies) and Measured 
!----------------------------------------------------

!-------------------------------------------------------------
! Graphics File:  Markers of the Bodies
!-------------------------------------------------------------

! Generic sphere for measured markers
Graphic=MeasuredMarker SPHERE x0=0.0  y0=0.0  z0=0.0  Radius=MarkerRadius;  Facets=40  Material=green

! MarkerSet_G10_RS graphics in Body Shank
Graphic=Shank_MarkerSet_G10_RS SPHERE x0=0.0301  y0=0.1213  z0=0.0157  Radius=MarkerRadius;  Facets=40  Material=blue
Object=Shank   Entity=Shank_MarkerSet_G10_RS   Motion=myMotion

! MarkerSet_G10_LM graphics in Body Shank
Graphic=Shank_MarkerSet_G10_LM SPHERE x0=-0.0165  y0=-0.0091  z0=0.0386  Radius=MarkerRadius;  Facets=40  Material=blue
Object=Shank   Entity=Shank_MarkerSet_G10_LM   Motion=myMotion

! MarkerSet_G10_MM graphics in Body Shank
Graphic=Shank_MarkerSet_G10_MM SPHERE x0=0.0165  y0=0.009  z0=-0.0386  Radius=MarkerRadius;  Facets=40  Material=blue
Object=Shank   Entity=Shank_MarkerSet_G10_MM   Motion=myMotion

! MarkerSet_G10_F2 graphics in Body Foot
Graphic=Foot_MarkerSet_G10_F2 SPHERE x0=0.0611  y0=0.0029  z0=0.0325  Radius=MarkerRadius;  Facets=40  Material=blue
Object=Foot   Entity=Foot_MarkerSet_G10_F2   Motion=myMotion

! MarkerSet_G10_F1 graphics in Body Foot
Graphic=Foot_MarkerSet_G10_F1 SPHERE x0=0.1307  y0=-0.0376  z0=-0.0059  Radius=MarkerRadius;  Facets=40  Material=blue
Object=Foot   Entity=Foot_MarkerSet_G10_F1   Motion=myMotion

! MarkerSet_G10_F3 graphics in Body Foot
Graphic=Foot_MarkerSet_G10_F3 SPHERE x0=0.0833  y0=-0.041  z0=0.0858  Radius=MarkerRadius;  Facets=40  Material=blue
Object=Foot   Entity=Foot_MarkerSet_G10_F3   Motion=myMotion

! MarkerSet_G10_RT graphics in Body Thigh
Graphic=Thigh_MarkerSet_G10_RT SPHERE x0=0.0982  y0=0.1247  z0=-0.006  Radius=MarkerRadius;  Facets=40  Material=blue
Object=Thigh   Entity=Thigh_MarkerSet_G10_RT   Motion=myMotion

! MarkerSet_G10_LFE graphics in Body Thigh
Graphic=Thigh_MarkerSet_G10_LFE SPHERE x0=0  y0=0.006  z0=0.055  Radius=MarkerRadius;  Facets=40  Material=blue
Object=Thigh   Entity=Thigh_MarkerSet_G10_LFE   Motion=myMotion

! MarkerSet_G10_MFE graphics in Body Thigh
Graphic=Thigh_MarkerSet_G10_MFE SPHERE x0=0  y0=-0.006  z0=-0.055  Radius=MarkerRadius;  Facets=40  Material=blue
Object=Thigh   Entity=Thigh_MarkerSet_G10_MFE   Motion=myMotion


!-------------------------------------------------------------
! Graphics File:  Markers measured with VICON
!-------------------------------------------------------------

! Spheres
Object=MeasuredMarker_MarkerSet_G10_RS  Entity=MeasuredMarker   Motion=myMotion
Object=MeasuredMarker_MarkerSet_G10_LM  Entity=MeasuredMarker   Motion=myMotion
Object=MeasuredMarker_MarkerSet_G10_MM  Entity=MeasuredMarker   Motion=myMotion
Object=MeasuredMarker_MarkerSet_G10_F1  Entity=MeasuredMarker   Motion=myMotion
Object=MeasuredMarker_MarkerSet_G10_F2  Entity=MeasuredMarker   Motion=myMotion
Object=MeasuredMarker_MarkerSet_G10_F3  Entity=MeasuredMarker   Motion=myMotion
Object=MeasuredMarker_MarkerSet_G10_RT  Entity=MeasuredMarker   Motion=myMotion
Object=MeasuredMarker_MarkerSet_G10_LFE  Entity=MeasuredMarker   Motion=myMotion
Object=MeasuredMarker_MarkerSet_G10_MFE  Entity=MeasuredMarker   Motion=myMotion
!Cylinder between each experimental skin marker and each body skin marker 
Graphic=Aux_Wireframe CYLINDER x0=0  y0=0  z0=0  x1=1  y1=0  z1=0  Radius=AuxWireRadius;  Facets=40  Material=red
Object=Aux_MarkerSet_G10_RS   Entity=Aux_Wireframe   Motion=myMotion

Object=Aux_MarkerSet_G10_LM   Entity=Aux_Wireframe   Motion=myMotion

Object=Aux_MarkerSet_G10_MM   Entity=Aux_Wireframe   Motion=myMotion

Object=Aux_MarkerSet_G10_F1   Entity=Aux_Wireframe   Motion=myMotion

Object=Aux_MarkerSet_G10_F2   Entity=Aux_Wireframe   Motion=myMotion

Object=Aux_MarkerSet_G10_F3   Entity=Aux_Wireframe   Motion=myMotion

Object=Aux_MarkerSet_G10_RT   Entity=Aux_Wireframe   Motion=myMotion

Object=Aux_MarkerSet_G10_LFE   Entity=Aux_Wireframe   Motion=myMotion

Object=Aux_MarkerSet_G10_MFE   Entity=Aux_Wireframe   Motion=myMotion

