classdef JOINT < handle
    %JOINT Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        Name    % Name of Joint        char
        Seg1    % Name of Segment 1    char
        Seg2    % Name of Segment 2    char
        Type    % Type of Joint        char
        F       % The force in the joint in child JCS       double[NFrames x 3]
        M       % The moment in the joint in child JCS      double[NFrames x 3]
        FGlob   % The force in the joint in Glob CS         double[NFrames x 3]
        MGlob   % The moment in the joint in Glob CS        double[NFrames x 3]
        Sensor  % Handle to the sensor           handle SENSOR
        Angles  %                                double[NFrames x 3]
        Anglesdot %                              double[NFrames x 3]
        Angles2dot%                              double[NFrames x 3]
    end
    
    methods
         function J = JOINT(Name,Type,Seg1,Seg2)
            J.Name = Name;
            J.Seg1 = Seg1;
            J.Seg2 = Seg2;
            J.Type = Type;
         end
         function JointElement = writeSubjectParsIntermed(J,docNode)
             JointElement = docNode.createElement('Joint');
             NAttribute = docNode.createAttribute('Name');
             TAttribute = docNode.createAttribute('Type');
             S1Attribute = docNode.createAttribute('Seg1');
             S2Attribute = docNode.createAttribute('Seg2');
             C1Attribute = docNode.createAttribute('Centre1');
             C2Attribute = docNode.createAttribute('Centre2');
             NAttribute.setValue(J.Name);
             TAttribute.setValue(J.Type);
             S1Attribute.setValue(J.Seg1);
             S2Attribute.setValue(J.Seg2);
             C1Attribute.setValue(J.Point1.Point.Name);
             C2Attribute.setValue(J.Point2.Point.Name);
             JointElement.setAttributeNode(NAttribute);
             JointElement.setAttributeNode(TAttribute);
             JointElement.setAttributeNode(S1Attribute);
             JointElement.setAttributeNode(S2Attribute);
             JointElement.setAttributeNode(C1Attribute);
             JointElement.setAttributeNode(C2Attribute);             
         end
    end
    
end

