%% PROJECT 2 VELOCITY ESTIMATION
tic
close all;
clear all;
clc;
addpath('../data')

%Change this for both dataset 1 and dataset 4. Do not use dataset 9.
datasetNum = 4;


[sampledData, sampledVicon, sampledTime] = init(datasetNum);
%% INITIALIZE CAMERA MATRIX AND OTHER NEEDED INFORMATION
Camera_Matrix=[311.0520,0,201.8724;0,311.3885,113.6210;0,0,1];
for n = 2:length(sampledData)
    %% Initalize Loop load images
    prev_images=sampledData(n-1).img;
    curr_images=sampledData(n).img;
    %% Detect good points
    prev=detectFASTFeatures(prev_images);
    %[prevFeatures, prevPoints] = extractFeatures(prev_images, prev);
    prevPoints=prev.selectStrongest(150);
    %% Initalize the tracker to the last frame.
    tracker=vision.PointTracker('MaxBidirectionalError',1);
    initialize(tracker,prevPoints.Location,prev_images);
    %% Find the location of the next points;
    [currPoints,currPointsVerify]=tracker(curr_images);
    %% Calculate velocity
    % Use a for loop
    CoordPrevpoints=[];
    CoordCurrpoints=[];
    VELOCITY=[];
    U=[];
    V=[];
    
    for m=1:length(prevPoints)
        dt=sampledData(n).t-sampledData(n-1).t;
        NormalisedPrevxyz=[prevPoints.Location(m,1);prevPoints.Location(m,2);1];
        FCPrev=pinv(Camera_Matrix)*NormalisedPrevxyz;
        CoordPrevpoints(m,1)=FCPrev(1,1);
        CoordPrevpoints(m,2)=FCPrev(2,1);
        NormalisedCurrxyz=[currPoints(m,1);currPoints(m,2);1];
        FCCurr=pinv(Camera_Matrix)*NormalisedCurrxyz;
        CoordCurrpoints(m,1)=FCCurr(1,1);
        CoordCurrpoints(m,2)=FCCurr(2,1);
        diffx=CoordCurrpoints(m,1)-CoordPrevpoints(m,1);
        diffy=CoordCurrpoints(m,2)-CoordPrevpoints(m,2);
        U(m)=diffx/dt;
        V(m)=diffy/dt;
        VELOCITY=[VELOCITY;U(m);V(m)];
    end
    %% Calculate Height
    [position, orientation, R_c2w] = estimatePose(sampledData, n);
    rbwn=eul2rotm(orientation);
    Hbwn=[rbwn,transpose(position);0,0,0,1];
    rcbn=[1.414/2,-1.414/2,0;-1.414/2,-1.414/2,0;0,0,-1];
    tcbn=[0.04;0.0;-0.03];
    Hcbn=[rcbn,tcbn;0,0,0,1];
    Hcwn=pinv(Hbwn*Hcbn);
    Z=Hcwn(3,4);
    functionxyz=[];
    for l=1:length(currPoints)
        tempfunction=[-1/Z,0,CoordCurrpoints(l,1)/Z,CoordCurrpoints(l,1)*CoordCurrpoints(l,2),-(1+CoordCurrpoints(l,1)^2),CoordCurrpoints(l,2);
        0,-1/Z,CoordCurrpoints(l,2)/Z,(1+CoordCurrpoints(l,2)^2),-CoordCurrpoints(l,1)*CoordCurrpoints(l,2),-CoordCurrpoints(l,1)];
        functionxyz=[functionxyz;tempfunction];
    end  
    %% RANSAC    
    % Write your own RANSAC implementation in the file velocityRANSAC
    optPo=[CoordCurrpoints(:,1),CoordCurrpoints(:,2)];
    e=0.8;
    VR = velocityRANSAC(VELOCITY,optPo,Z,R_c2w,e);
    
    %% Thereshold outputs into a range.
    % Not necessary
    
    %% Fix the linear velocity
    % Change the frame of the computed velocity to world frame
    velocity=pinv(functionxyz)*VELOCITY;
    Hcwn=pinv(Hcwn);
    Vi=[Hcwn(1:3,1:3),-Hcwn(1:3,1:3)*skew(Hcwn(1:3,4));zeros(3),Hcwn(1:3,1:3)]*velocity;
    estimatedV(:,n)=Vi;
    %% ADD SOME LOW PASS FILTER CODE
    % Not neceessary but recommended 
    estimatedV(:,n) =Vi;
    %% STORE THE COMPUTED VELOCITY IN THE VARIABLE estimatedV AS BELOW
   % estimatedV(:,n) = Vel; % Feel free to change the variable Vel to anything that you used.
    % Structure of the Vector Vel should be as follows:
    % Vel(1) = Linear Velocity in X
    % Vel(2) = Linear Velocity in Y
    % Vel(3) = Linear Velocity in Z
    % Vel(4) = Angular Velocity in X
    % Vel(5) = Angular Velocity in Y
    % Vel(6) = Angular Velocity in Z
end
    time_vector= zeros(length(sampledData), 1);
        for tv=1:length(time_vector)
            time_vector(tv,1) = sampledData(1,tv).t;
        end
    time_vector = sgolayfilt(time_vector, 1, 101);

    vx = transpose(estimatedV(1,:));
    vy = transpose(estimatedV(2,:));
    vz = transpose(estimatedV(3,:));

    wx = transpose(estimatedV(4,:));
    wy = transpose(estimatedV(5,:));
    wz = transpose(estimatedV(6,:));


    estimatedV(1,:) = sgolayfilt(vx, 1, 17);
    estimatedV(2,:) = sgolayfilt(vy, 1, 17);
    estimatedV(3,:) = sgolayfilt(vz, 1, 7);
    estimatedV(4,:) = sgolayfilt(wx, 1, 5);
    estimatedV(5,:) = sgolayfilt(wy, 1, 5);
    estimatedV(6,:) = sgolayfilt(wz, 1, 5);
    
plotData(estimatedV, sampledData, sampledVicon, sampledTime, datasetNum)
toc