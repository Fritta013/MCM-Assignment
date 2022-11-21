%% Exercises Modelling Part 1
% Rotation matrices, Equivalent angle-axis representations, Quaternions
addpath('include') %%DO NOT CHANGE STUFF INSIDE THIS PATH

%% Exercise 1
% Implementing the Rodrigues Formula, 
% taking in input the geometric unit vector v and the rotation angle theta
% and returning the orientation matrix.
% testing for following cases:

%% 1.1.
%Test first for initial frame, no rotation 
theta = 0; 
v=[1, 0, 0]; 
aRb = ComputeAngleAxis(theta, v)
disp('aRb ex 1.1:');disp(aRb);
plotRotation(theta,v,aRb);
disp('theta ex 1.2:');disp(theta);
disp('v ex 1.2:');disp(v);

%Rotation around x axis for 30 deg
theta = deg2rad(30); 
v=[1, 0, 0]; 
aRb = ComputeAngleAxis(theta, v)
disp('aRb ex 1.1:');disp(aRb);
plotRotation(theta,v,aRb);
disp('theta ex 1.2:');disp(theta);
disp('v ex 1.2:');disp(v);

%% 1.3.
%Rotation around y axis for angle pi/4
theta = pi/4; 
v=[0, 1, 0]; 
aRb = ComputeAngleAxis(theta, v)
disp('aRb ex 1.1:');disp(aRb);
plotRotation(theta,v,aRb);
disp('theta ex 1.2:');disp(theta);
disp('v ex 1.2:');disp(v);

%% 1.4
%Rotation around z axis for angle pi/2
theta = pi/2; 
v=[0, 0, 1]; 
aRb = ComputeAngleAxis(theta, v)
disp('aRb ex 1.1:');disp(aRb);
plotRotation(theta,v,aRb);
disp('theta ex 1.2:');disp(theta);
disp('v ex 1.2:');disp(v);

%% 1.5.

theta = 0.2449; 
v=[0.408, 0.816, -0.408]; 
aRb = ComputeAngleAxis(theta, v)
disp('aRb ex 1.1:');disp(aRb);
plotRotation(theta,v,aRb);
disp('theta ex 1.2:');disp(theta);
disp('v ex 1.2:');disp(v);


%% 1.6.
rho=[0, pi/2, 0];
theta = norm(rho) %to obtain angle of rotation
v=rho/theta %to obtain axis of rotation
aRb = ComputeAngleAxis(theta, v)
disp('aRb ex 1.1:');disp(aRb);
plotRotation(theta,v,aRb);
disp('theta ex 1.2:');disp(theta);
disp('v ex 1.2:');disp(v);

%% 1.7.
rho= [0.4, -0.3, -0.3]; 
theta = norm(rho) %to obtain angle of rotation
v=rho/theta %to obtain axis of rotation
aRb = ComputeAngleAxis(theta, v)
disp('aRb ex 1.1:');disp(aRb);
plotRotation(theta,v,aRb);
disp('theta ex 1.2:');disp(theta);
disp('v ex 1.2:');disp(v);

%% 1.8.
 
rho=  [-pi/4, -pi/3, pi/8];
theta = norm(rho) %to obtain angle of rotation
v=rho/theta %to obtain axis of rotation
aRb = ComputeAngleAxis(theta, v)
disp('aRb ex 1.1:');disp(aRb);
plotRotation(theta,v,aRb);
disp('theta ex 1.2:');disp(theta);
disp('v ex 1.2:');disp(v);

%% Exercise 2

%% 2.1 

%from figure 1, the rotation axis is z and the rotation angle is pi/2, then
v=[0,0,1]
theta=pi/2
% Compute the rotation matrix between frame <a> and <b>
%From the orientation matrix developed in exercice 1, we find aRb
aRb= ComputeAngleAxis(theta, v)

% 2.2. Solve the Inverse Equivalent Angle-Axis Problem for the orientation matrix aRb. 
% Compute the inverse equivalent angle-axis repr. of aRb 
    [theta, v] = ComputeInverseAngleAxis(aRb)
    % Plot Results
    disp('aRb ex 2.1:');disp(aRb);
    plotRotation(theta,v,aRb);
    disp('theta ex 2.2:');disp(theta);
    disp('v ex 2.2:');disp(v); 


 %% 2.3
    wRc = [0.835959, -0.283542, -0.469869; 
        0.271321, 0.957764, -0.0952472;
        0.47703, -0.0478627, 0.877583]
%compute the R of c with respect to W then extract the angles
%finding the Inverse Equivalent Angle-Axis Problem for wRc, to find the
%rotation axis and angle
   [theta, v] = ComputeInverseAngleAxis(wRc);
%According to the obtained v and theta the rotation axis is and 
   thetab = pi/2; 
   vb=[ 0 , 0, 1]; 
   wRb = ComputeAngleAxis(thetab, vb);
    % Compute the rotation matrix between frame <c> and <b>
    cRb=transpose(wRc)*wRb
    % Compute inverse equivalent angle-axis repr. of cRb
    [theta, v] = ComputeInverseAngleAxis(cRb)
    % Plot Results
    plotRotation(theta,v,cRb);
    disp('theta ex 2.3:');disp(theta);
    disp('v ex 2.3:');disp(v); 

%% Exercise 3
%% 3.1
% a
        alpha=((30*pi)/180); %  angle of rotation 
    
        %rotation matrix from <w> to frame <b> by rotating around z-axes
        wRb_z = [cos(alpha), -sin(alpha), 0; sin(alpha) , cos(alpha), 0; 0, 0, 1]
% b
        beta=(45*pi)/180; %  angle of rotation 
    
        %rotation matrix from <w> to frame <b> by rotating around y-axes
        wRb_y = [cos(beta) , 0 , sin(beta) ; 0 , 1 , 0 ; -sin(beta) , 0 , cos(beta)] 
%c
       theta=(15*pi)/180; %  angle of rotation 

        %rotation matrix from <w> to frame <b> by rotating around x-axes
        wRb_x = [ 1 , 0 , 0 ; 0 , cos(theta) , -sin(theta) ; 0 , sin(theta) , cos(theta)] 
        disp('es 3.1:');disp(wRb_z);disp(wRb_y);disp(wRb_x);
    

    % a, rotation around z axis then 
        v=[0,0,1];
        %[alpha, v] = ComputeInverseAngleAxis(wRb_z);
        wRb_z = ComputeAngleAxis(alpha,v)
        % Plot Results
        plotRotation(alpha,v,wRb_z);
        disp('theta ex 3.2.a:');disp(alpha);
        disp('v ex 3.2.a:');disp(v); 

    % b, rotation around y axis then 
        v=[0,1,0];
        %[beta, v] = ComputeInverseAngleAxis(wRb_y);
        wRb_y = ComputeAngleAxis(beta,v);
        % Plot Results
        plotRotation(beta,v,wRb_y);
        disp('theta ex 3.2.b:');disp(beta);
        disp('v ex 3.2.b:');disp(v); 

    % c, , rotation around x axis then 
    v=[1,0,0];
        %[theta, v] = ComputeInverseAngleAxis(wRb_x);
        wRb_x = ComputeAngleAxis(theta,v);
        % Plot Results
        plotRotation(theta,v,wRb_x);
        disp('theta ex 3.2.c:');disp(theta);
        disp('v ex 3.2.c:');disp(v); 
%% 3.2 
    % Compute the rotation matrix corresponding to the z-y-x representation
    %Rzyx=Rot(z,alpha)*Rot(y,beta)*Rot(x,theta), then
    Rzyx=wRb_z*wRb_y*wRb_x
    % Compute equivalent angle-axis represenation
    [theta, v] = ComputeInverseAngleAxis(Rzyx)
    % Plot Results
    plotRotation(theta,v,Rzyx);
    disp('theta ex 3.3:');disp(theta);
    disp('v ex 3.3:');disp(v);
%% 3.3
    % Compute the rotation matrix corresponding to the z-x-z representation;
    %Rzxz=Rot(z,alpha)*Rot(x,theta)*Rot(z,alpha), then
    Rzxz = wRb_z*wRb_x*wRb_z
	% Compute inverse equivalent angle-axis repr.
    [theta, v] = ComputeInverseAngleAxis(Rzxz)
    % Plot Results
    plotRotation(theta,v,Rzxz);
    disp('theta ex 3.4:');disp(theta);
    disp('v ex 3.4:');disp(v);
%% Exercise 4
%% 4.1 
% we have q = 0.8924 +  0.23912i +  0.36964j + 0.099046k
    %%%%%%%%%%%%
    a = 0.8924; % real part % check notazione usata in classe
    b = 0.23912; %imaginary part 
    c = 0.36964; %imaginary part 
    d = 0.099046; %imaginary part 
    %%%%%%%%%%%%%
    % Compute the rotation matrix associated with the given quaternion
    rotMatrix = quatToRot(a,b,c,d) %return corresponding rotation matrix
    disp('rot matrix es 4.1');disp(rotMatrix)

    %% 4.1
    % solve using matlab functions quaternion(), rotmat(),
    rotMatrixMatlab =  quaternion(a,b,c,d) %create a quaternion array  
    rotationMatrix = rotmat(rotMatrixMatlab, 'point') %return corresponding rotation matrix 
    % same results as before

    [theta, v] = ComputeInverseAngleAxis(rotMatrix) %Inverse Equivalent Angle-Axis 
    % Plot Results
    plotRotation(theta,v,rotMatrix);
    disp('rot matrix es 4.2');disp(rotMatrix)