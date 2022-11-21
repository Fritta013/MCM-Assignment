function [rot_matrix] = quatToRot(q0,q1,q2,q3)
% quatToRot convert a quaternion into a rotation matrix
    %Convert a quaternion into a full three-dimensional rotation matrix.
    %Input
    %:param Q: A 4 element array representing the quaternion (q0,q1,q2,q3) 
    Q= [q0,q1,q2,q3]
    %Output
    %return: A 3x3 element matrix representing the full 3D rotation matrix. 
    
    %First row of the rotation matrix
    Rq1=[((2*((q0^2)+(q1^2)))-1) , (2*((q1*q2)-(q0*q3))) ,  (2*((q1*q3)+(q0*q2))) ]
    
    %Second row of the rotation matrix
    Rq2=[(2*((q1*q2)+(q0*q3))) , ((2*((q0^2)+(q2^2)))-1) ,  (2*((q2*q3)-(q0*q1))) ]

    %Third row of the rotation matrix
    Rq3=[(2*((q1*q3)-(q0*q2))) ,  (2*((q2*q3)+(q0*q1)))  ,  ((2*((q0^2)+(q3^2)))-1) ]
    %3x3 rotation matrix
    rot_matrix = [ Rq1 ; Rq2 ; Rq3]
end