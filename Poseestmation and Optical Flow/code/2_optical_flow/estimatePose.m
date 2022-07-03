function [position, orientation, R_c2w] = estimatePose(data,t)
%% CHANGE THE NAME OF THE FUNCTION TO estimatePose
% Please note that the coordinates for each corner of each AprilTag are
% defined in the world frame, as per the information provided in the
% handout. Ideally a call to the function getCorner with ids of all the
% detected AprilTags should be made. This function should return the X and
% Y coordinate of each corner, or each corner and the centre, of all the
% detected AprilTags in the image. You can implement that anyway you want
% as long as the correct output is received. A call to that function
% should made from this function.
id=data(t).id;
res=getCorner(id);
pimagezero=data(t).p0;
pimageone=data(t).p1;
pimagetwo=data(t).p2;
pimagethree=data(t).p3;
pimagefour=data(t).p4;
px=[pimagezero(1,:),pimageone(1,:),pimagetwo(1,:),pimagethree(1,:),pimagefour(1,:)]';
py=[pimagezero(2,:),pimageone(2,:),pimagetwo(2,:),pimagethree(2,:),pimagefour(2,:)]';
A=[];
for f=1:length(id)
    a=[res(f,1),res(f,2),1,0,0,0,-pimagezero(1,f)*res(f,1),-pimagezero(1,f)*res(f,2),-pimagezero(1,f);
        0,0,0,res(f,1),res(f,2),1,-pimagezero(2,f)*res(f,1),-pimagezero(2,f)*res(f,2),-pimagezero(2,f);
        
        res(length(id)+f,1),res(length(id)+f,2),1,0,0,0,-pimageone(1,f)*res(length(id)+f,1),-pimageone(1,f)*res(length(id)+f,2),-pimageone(1,f);
        0,0,0,res(length(id)+f,1),res(length(id)+f,2),1,-pimageone(2,f)* res(length(id)+f,1),-pimageone(2,f)*res(length(id)+f,2),-pimageone(2,f);
        
        res(2*length(id)+f,1),res(2*length(id)+f,2),1,0,0,0,-pimagetwo(1,f)*res(2*length(id)+f,1),-pimagetwo(1,f)*res(2*length(id)+f,2),-pimagetwo(1,f);
        0,0,0,res(2*length(id)+f,1),res(2*length(id)+f,2),1,-pimagetwo(2,f)*res(2*length(id)+f,1),-pimagetwo(2,f)*res(2*length(id)+f,2),-pimagetwo(2,f);
        
        res(3*length(id)+f,1),res(3*length(id)+f,2),1,0,0,0,-pimagethree(1,f)*res(3*length(id)+f,1),-pimagethree(1,f)*res(3*length(id)+f,2),-pimagethree(1,f);
        0,0,0,res(3*length(id)+f,1),res(3*length(id)+f,2),1,-pimagethree(2,f)*res(3*length(id)+f,1),-pimagethree(2,f)*res(3*length(id)+f,2),-pimagethree(2,f);
        
        res(4*length(id)+f,1),res(4*length(id)+f,2),1,0,0,0,-pimagefour(1,f)*res(4*length(id)+f,1),-pimagefour(1,f)*res(4*length(id)+f,2),-pimagefour(1,f);
        0,0,0,res(4*length(id)+f,1),res(4*length(id)+f,2),1,-pimagefour(2,f)*res(4*length(id)+f,1),-pimagethree(2,f)*res(4*length(id)+f,2),-pimagefour(2,f)];
    A=[A;a];
end
[u,s,v]=svd(A);
h=reshape(v(:,9),[3,3]);
H=transpose(h)/h(3,3);
Camera_Matrix=[311.0520,0,201.8724;0,311.3885,113.6210;0,0,1];
rt=pinv(Camera_Matrix)*H;
rone=rt(:,1);
norm(rone)
rtwo=rt(:,2);
b=rt(:,3);
rcot=cross(rone,rtwo);
[uo,so,vo]=svd([rone,rtwo,rcot]);
rwc=uo*[1,0,0;0,1,0;0,0,det(uo*transpose(vo))]*transpose(vo);
twc=b/norm(rone);
Hwc=[rwc,twc;0,0,0,1];
rcb=[1.414/2,-1.414/2,0;-1.414/2,-1.414/2,0;0,0,-1];
tcb=[0.04;0.0;-0.03];
Hcb=[rcb,tcb;0,0,0,1];
hbw=Hcb*Hwc;
Hbw=pinv(hbw);
rbw=Hbw(1:3,1:3);
tbw=Hbw(1:3,4);
oa=rotm2eul(rbw,'ZYX');
R_c2w=pinv(rwc);
%T_c2w=-pinv(rwc)*twc;
position=transpose(tbw);
orientation=oa;

    %% Input Parameter Defination
    % data = the entire data loaded in the current dataset
    % t = index of the current data in the dataset
    
    %% Output Parameter Defination
    
    % position = translation vector representing the position of the
    % drone(body) in the world frame in the current time, in the order XYZ
    
    % orientation = euler angles representing the orientation of the
    % drone(body) in the world frame in the current time, in the order XYZ
    
    %R_c2w = Rotation which defines camera to world frame
end