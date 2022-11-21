function R = ComputeAngleAxis(theta,v)

%extracting x, y, z from v 
vx1=v(1);
vy2=v(2);
vz3=v(3);
k = [ 0 -vz3 vy2; vz3 0 -vx1 ; -vy2 vx1 0];
%calculate rotation matrix based on Rodrigues formula
I = eye(3); %idendity matrix 
R= [I+(sin(theta)*k)+((1-cos(theta))*(k*k))];

end 