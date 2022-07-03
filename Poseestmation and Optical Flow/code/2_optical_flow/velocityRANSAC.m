function [Vel] = velocityRANSAC(optV,optPos,Z,R_c2w,e)
%% CHANGE THE NAME OF THE FUNCTION TO velocityRANSAC
    %% Input Parameter Description
    % optV = The optical Flow
    % optPos = Position of the features in the camera frame 
    % Z = Height of the drone
    % R_c2w = Rotation defining camera to world frame
    % e = RANSAC hyper parameter    
    p_success=0.99;
    Homography_M=4;
    e=0.8;
    max_no_iterations=log(1-p_success)/log(1-e^Homography_M);
    count=0;
    for iteration=1:max_no_iterations
        liners=0;
        points=randperm(length(optPos),3);
        point_one=optPos(points(1,1),:);
        point_two=optPos(points(1,2),:);
        point_three=optPos(points(1,3),:);
        h=[-1/Z,0,point_one(1,1)/Z,point_one(1,1)*point_one(1,2),-(1+point_one(1,1)^2),point_one(1,2);
            0,-1/Z,point_one(1,2)/Z,(1+point_one(1,2)^2),-point_one(1,1)*point_one(1,2),-point_one(1,1);
            -1/Z,0,point_two(1,1)/Z,point_two(1,1)*point_two(1,2),-(1+point_two(1,1)^2),point_two(1,2);
            0,-1/Z,point_two(1,2)/Z,(1+point_two(1,2)^2),-point_two(1,1)*point_two(1,2),-point_two(1,1);
            -1/Z,0,point_three(1,1)/Z,point_three(1,1)*point_three(1,2),-(1+point_three(1,1)^2),point_three(1,2);
            0,-1/Z,point_three(1,2)/Z,(1+point_three(1,2)^2),-point_three(1,1)*point_three(1,2),-point_three(1,1)];
        Voptimal=[optV(2*points(1,1)-1);optV(2*points(1,1));optV(2*points(1,2)-1);optV(2*points(1,2));optV(2*points(1,3)-1);optV(2*points(1,3))];
        Vj=pinv(h)*Voptimal;
        for x=1:length(optPos)
            tempH=[-1/Z,0,optPos(x,1)/Z,optPos(x,1)*optPos(x,2),-(1+optPos(x,1)^2),optPos(x,2);
                0,-1/Z,optPos(x,2)/Z,(1+optPos(x,2)^2),-optPos(x,1)*optPos(x,2),-optPos(x,1)];
            tempP=[optV(2*x-1);optV(2*x)];
            change=norm(tempH*Vj-tempP)^2;
            if(change<=1e-3)
                liners=liners+1;
            end
        end
        if(liners>=count)
             count = liners;
             Vel=Vj;
        end
   end
    
    %% Output Parameter Description
    % Vel = Linear velocity and angualr velocity vector
    
end