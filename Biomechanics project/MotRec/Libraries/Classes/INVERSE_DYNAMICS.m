classdef INVERSE_DYNAMICS< handle
    %INVERSE_DYNAMIC Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        
    end
    
    methods
        function ID = INVERSE_DYNAMICS()
        end
        function [Child_w,Child_wdot,Child_vdotOrBCS,Child_F_CoM,Child_M_CoM]= calcInertiaForces(ID,...
                Child_R_Parent,Parent_w,Parent_wdot,Parent2Child_Anglesdot,Parent2Child_Angles2dot,Parent_PtJnt2Child,Parent_vdotOrBCS,Child_PtOrBCS,Child_PtCoM,Child_m,ChildCoM_I) %#ok<MANU>
            Child_w    = Child_R_Parent * Parent_w + Parent2Child_Anglesdot;
            Child_wdot = Child_R_Parent * Parent_wdot + cross(Child_R_Parent * Parent_w,Parent2Child_Anglesdot) + Parent2Child_Angles2dot;
            Child_vdotJnt2Par = Child_R_Parent *(cross(Parent_wdot,Parent_PtJnt2Child)+cross(Parent_w,cross(Parent_w,Parent_PtJnt2Child))+ Parent_vdotOrBCS);
            Child_vdotOrBCS = Child_vdotJnt2Par + cross(Child_wdot,Child_PtOrBCS)+cross(Child_w,cross(Child_w,Child_PtOrBCS));
            Child_vdotCoM = Child_vdotOrBCS + cross(Child_wdot,Child_PtCoM)+cross(Child_w,cross(Child_w,Child_PtCoM));
            Child_F_CoM = Child_m*Child_vdotCoM;
            Child_M_CoM = ChildCoM_I * Child_wdot + cross(Child_w,ChildCoM_I*Child_w);
        end
%         function [Segment_F_CoM,Segment_M_CoM] = calcInertiaFoces(ID,Segment_w,Segment_wdot,Segment_PtCoM,Segment_vdotOr,...
%                                                 Segment_m,SegmentCoM_I) %#ok<MANU>
%             Segment_vdotCoM = cross(Segment_wdot,Segment_PtCoM)+cross(Segment_w,cross(Segment_w,Segment_PtCoM))+ Segment_vdotOr;
%             Segment_F_CoM = Segment_m*Segment_vdotCoM;
%             Segment_M_CoM = SegmentCoM_I * Segment_wdot + cross(Segment_w,SegmentCoM_I*Segment_w);
%             
%         end
%         function [Parent_F_Jnt2Parent,Parent_M_Jnt2Parent] = calcJointForces(ID,Parent_F_CoM,Parent_PtCoM,Parent_M_CoM,...
%                 Parent_F_Ext,Parent_PtExt,Parent_M_Ext,Paren_R_Child,Child_F_Jnt2Parent,Parent_PtJnt2Child,Child_M_Jnt2Parent)  %#ok<MANU>
%             Parent_F_Jnt2Parent = Parent_F_CoM + Parent_F_Ext + Paren_R_Child * Child_F_Jnt2Parent;
%             Parent_M_Jnt2Parent = cross(Parent_PtCoM,Parent_F_CoM)+ Parent_M_CoM + cross(Parent_PtExt,Parent_F_Ext)+ ...
%                 Parent_M_Ext + cross(Parent_PtJnt2Child,Paren_R_Child*Child_F_Jnt2Parent)+ Paren_R_Child * Child_M_Jnt2Parent;
%                 
%         end
        function [Seg_F_Jnt2Parent,Seg_M_Jnt2Parent] = calcJointForces(ID,Seg_F_CoM,Seg_F_Ext,Seg_F_Jnt2Child,...
                Seg_M_Jnt2Child,Seg_M_Ext,Seg_M_CoM,Seg_Fd_Jnt2Child,Seg_Fd_Ext,Seg_Fd_CoM,Seg_Pos_Jnt2Parent)  %#ok<MANU>
            Seg_F_Jnt2Parent = Seg_F_CoM - Seg_F_Ext - Seg_F_Jnt2Child;
            Seg_M_Jnt2Parent = Seg_M_CoM -Seg_M_Jnt2Child - Seg_M_Ext - Seg_Fd_Jnt2Child - Seg_Fd_Ext + Seg_Fd_CoM...
                               - cross(Seg_Pos_Jnt2Parent,Seg_F_Jnt2Parent);                
        end
    end
    
end

