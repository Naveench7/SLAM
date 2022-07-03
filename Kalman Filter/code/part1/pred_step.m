 function [covarEst,uEst] = pred_step(uPrev,covarPrev,angVel,acc,dt)
%covarPrev and uPrev are the previous mean and covariance respectively
%angVel is the angular velocity
%acc is the accelerationsy
%dt is the sampling time
gx=0;
gy=0;
gz=9.81;
vx=uPrev(7);
vy=uPrev(8);
vz=uPrev(9);
orix=uPrev(4);
oriy=uPrev(5);
oriz=uPrev(6);
bgx=uPrev(10);
bgy=uPrev(11);
bgz=uPrev(12);
bax=uPrev(13);
bay=uPrev(14);
baz=uPrev(15);
omegamx=angVel(1);
omegamy=angVel(2);
omegamz=angVel(3);
amx=acc(1);
amy=acc(2);
amz=acc(3);
ngx=0;
ngy=0;
ngz=0;
nax=0;
nay=0;
naz=0;
nbgx=1.2;
nbgy=1.5;
nbgz=1.5;
nbax=1.5;
nbay=1.5;
nbaz=1.5;
V=[vx;vy;vz];
g=[gx;gy;gz];
G=[0 -sin(oriz) cos(oriz)*cos(oriy);0 cos(oriz) sin(oriz)*cos(oriy);1 0 -sin(oriy)];
Gi=inv(G);
Omegam=[omegamx;omegamy;omegamz];
Bg=[bgx;bgy;bgz];
Ng=[ngx;ngy;ngz];
Rz=[cos(oriz) -sin(oriz) 0;sin(oriz) cos(oriz) 0;0 0 1];
Ry=[cos(oriy) 0 sin(oriy);0 1 0;-sin(oriy) 0 cos(oriy)];
Rx=[1 0 0;0 cos(orix) -sin(orix);0 sin(orix) cos(orix)];
R=Rz*Ry*Rx;
Am=[amx;amy;amz];
Ba=[bax;bay;baz];
Na=[nax;nay;naz];
Nbg=[nbgx;nbgy;nbgz];
Nba=[nbax;nbay;nbaz];
x_dot=[V;Gi*R*(Omegam-Bg-Ng);g+(R*(Am-Ba-Na));Nbg;Nba];
At=[0,0,0,0,0,0,1,0,0,0,0,0,0,0,0;0,0,0,0,0,0,0,1,0,0,0,0,0,0,0;0,0,0,0,0,0,0,0,1,0,0,0,0,0,0;0,0,0,(cos(oriy)*sin(orix) - (cos(oriz)*sin(oriy)*(cos(orix)*sin(oriz) - cos(oriz)*sin(orix)*sin(oriy)))/(cos(oriy)*cos(oriz)^2 + cos(oriy)*sin(oriz)^2) + (sin(oriy)*sin(oriz)*(cos(orix)*cos(oriz) + sin(orix)*sin(oriy)*sin(oriz)))/(cos(oriy)*cos(oriz)^2 + cos(oriy)*sin(oriz)^2))*(bgz + ngz - omegamz) - (cos(orix)*cos(oriy) + (cos(oriz)*sin(oriy)*(sin(orix)*sin(oriz) + cos(orix)*cos(oriz)*sin(oriy)))/(cos(oriy)*cos(oriz)^2 + cos(oriy)*sin(oriz)^2) - (sin(oriy)*sin(oriz)*(cos(oriz)*sin(orix) - cos(orix)*sin(oriy)*sin(oriz)))/(cos(oriy)*cos(oriz)^2 + cos(oriy)*sin(oriz)^2))*(bgy + ngy - omegamy), - (bgx + ngx - omegamx)*((cos(oriy)^2*cos(oriz)^2)/(cos(oriy)*cos(oriz)^2 + cos(oriy)*sin(oriz)^2) - cos(oriy) + (cos(oriy)^2*sin(oriz)^2)/(cos(oriy)*cos(oriz)^2 + cos(oriy)*sin(oriz)^2) - (cos(oriz)^2*sin(oriy)^2)/(cos(oriy)*cos(oriz)^2 + cos(oriy)*sin(oriz)^2) - (sin(oriy)^2*sin(oriz)^2)/(cos(oriy)*cos(oriz)^2 + cos(oriy)*sin(oriz)^2) + (cos(oriy)*cos(oriz)^2*sin(oriy)*(sin(oriy)*cos(oriz)^2 + sin(oriy)*sin(oriz)^2))/(cos(oriy)*cos(oriz)^2 + cos(oriy)*sin(oriz)^2)^2 + (cos(oriy)*sin(oriy)*sin(oriz)^2*(sin(oriy)*cos(oriz)^2 + sin(oriy)*sin(oriz)^2))/(cos(oriy)*cos(oriz)^2 + cos(oriy)*sin(oriz)^2)^2) - (bgz + ngz - omegamz)*((cos(oriy)*cos(oriz)*(sin(orix)*sin(oriz) + cos(orix)*cos(oriz)*sin(oriy)))/(cos(oriy)*cos(oriz)^2 + cos(oriy)*sin(oriz)^2) - cos(orix)*sin(oriy) - (cos(oriy)*sin(oriz)*(cos(oriz)*sin(orix) - cos(orix)*sin(oriy)*sin(oriz)))/(cos(oriy)*cos(oriz)^2 + cos(oriy)*sin(oriz)^2) + (cos(oriz)*sin(oriy)*(sin(oriy)*cos(oriz)^2 + sin(oriy)*sin(oriz)^2)*(sin(orix)*sin(oriz) + cos(orix)*cos(oriz)*sin(oriy)))/(cos(oriy)*cos(oriz)^2 + cos(oriy)*sin(oriz)^2)^2 - (sin(oriy)*sin(oriz)*(sin(oriy)*cos(oriz)^2 + sin(oriy)*sin(oriz)^2)*(cos(oriz)*sin(orix) - cos(orix)*sin(oriy)*sin(oriz)))/(cos(oriy)*cos(oriz)^2 + cos(oriy)*sin(oriz)^2)^2 + (cos(orix)*cos(oriy)*cos(oriz)^2*sin(oriy))/(cos(oriy)*cos(oriz)^2 + cos(oriy)*sin(oriz)^2) + (cos(orix)*cos(oriy)*sin(oriy)*sin(oriz)^2)/(cos(oriy)*cos(oriz)^2 + cos(oriy)*sin(oriz)^2)) - (bgy + ngy - omegamy)*((cos(oriy)*sin(oriz)*(cos(orix)*cos(oriz) + sin(orix)*sin(oriy)*sin(oriz)))/(cos(oriy)*cos(oriz)^2 + cos(oriy)*sin(oriz)^2) - (cos(oriy)*cos(oriz)*(cos(orix)*sin(oriz) - cos(oriz)*sin(orix)*sin(oriy)))/(cos(oriy)*cos(oriz)^2 + cos(oriy)*sin(oriz)^2) - sin(orix)*sin(oriy) - (cos(oriz)*sin(oriy)*(sin(oriy)*cos(oriz)^2 + sin(oriy)*sin(oriz)^2)*(cos(orix)*sin(oriz) - cos(oriz)*sin(orix)*sin(oriy)))/(cos(oriy)*cos(oriz)^2 + cos(oriy)*sin(oriz)^2)^2 + (sin(oriy)*sin(oriz)*(sin(oriy)*cos(oriz)^2 + sin(oriy)*sin(oriz)^2)*(cos(orix)*cos(oriz) + sin(orix)*sin(oriy)*sin(oriz)))/(cos(oriy)*cos(oriz)^2 + cos(oriy)*sin(oriz)^2)^2 + (cos(oriy)*cos(oriz)^2*sin(orix)*sin(oriy))/(cos(oriy)*cos(oriz)^2 + cos(oriy)*sin(oriz)^2) + (cos(oriy)*sin(orix)*sin(oriy)*sin(oriz)^2)/(cos(oriy)*cos(oriz)^2 + cos(oriy)*sin(oriz)^2)),0,0,0,0,0,0,0,0,0,0;0,0,0,((cos(oriz)*(cos(oriz)*sin(orix) - cos(orix)*sin(oriy)*sin(oriz)))/(cos(oriz)^2 + sin(oriz)^2) + (sin(oriz)*(sin(orix)*sin(oriz) + cos(orix)*cos(oriz)*sin(oriy)))/(cos(oriz)^2 + sin(oriz)^2))*(bgy + ngy - omegamy) + ((cos(oriz)*(cos(orix)*cos(oriz) + sin(orix)*sin(oriy)*sin(oriz)))/(cos(oriz)^2 + sin(oriz)^2) + (sin(oriz)*(cos(orix)*sin(oriz) - cos(oriz)*sin(orix)*sin(oriy)))/(cos(oriz)^2 + sin(oriz)^2))*(bgz + ngz - omegamz),0,0,0,0,0,0,0,0,0,0,0;0,0,0,- ((cos(oriz)*(sin(orix)*sin(oriz) + cos(orix)*cos(oriz)*sin(oriy)))/(cos(oriy)*cos(oriz)^2 + cos(oriy)*sin(oriz)^2) - (sin(oriz)*(cos(oriz)*sin(orix) - cos(orix)*sin(oriy)*sin(oriz)))/(cos(oriy)*cos(oriz)^2 + cos(oriy)*sin(oriz)^2))*(bgy + ngy - omegamy) - ((cos(oriz)*(cos(orix)*sin(oriz) - cos(oriz)*sin(orix)*sin(oriy)))/(cos(oriy)*cos(oriz)^2 + cos(oriy)*sin(oriz)^2) - (sin(oriz)*(cos(orix)*cos(oriz) + sin(orix)*sin(oriy)*sin(oriz)))/(cos(oriy)*cos(oriz)^2 + cos(oriy)*sin(oriz)^2))*(bgz + ngz - omegamz),(bgx + ngx - omegamx)*((cos(oriz)^2*sin(oriy))/(cos(oriy)*cos(oriz)^2 + cos(oriy)*sin(oriz)^2) + (sin(oriy)*sin(oriz)^2)/(cos(oriy)*cos(oriz)^2 + cos(oriy)*sin(oriz)^2) - (cos(oriy)*cos(oriz)^2*(sin(oriy)*cos(oriz)^2 + sin(oriy)*sin(oriz)^2))/(cos(oriy)*cos(oriz)^2 + cos(oriy)*sin(oriz)^2)^2 - (cos(oriy)*sin(oriz)^2*(sin(oriy)*cos(oriz)^2 + sin(oriy)*sin(oriz)^2))/(cos(oriy)*cos(oriz)^2 + cos(oriy)*sin(oriz)^2)^2) - (bgz + ngz - omegamz)*((cos(oriz)*(sin(oriy)*cos(oriz)^2 + sin(oriy)*sin(oriz)^2)*(sin(orix)*sin(oriz) + cos(orix)*cos(oriz)*sin(oriy)))/(cos(oriy)*cos(oriz)^2 + cos(oriy)*sin(oriz)^2)^2 - (sin(oriz)*(sin(oriy)*cos(oriz)^2 + sin(oriy)*sin(oriz)^2)*(cos(oriz)*sin(orix) - cos(orix)*sin(oriy)*sin(oriz)))/(cos(oriy)*cos(oriz)^2 + cos(oriy)*sin(oriz)^2)^2 + (cos(orix)*cos(oriy)*cos(oriz)^2)/(cos(oriy)*cos(oriz)^2 + cos(oriy)*sin(oriz)^2) + (cos(orix)*cos(oriy)*sin(oriz)^2)/(cos(oriy)*cos(oriz)^2 + cos(oriy)*sin(oriz)^2)) - (bgy + ngy - omegamy)*((sin(oriz)*(sin(oriy)*cos(oriz)^2 + sin(oriy)*sin(oriz)^2)*(cos(orix)*cos(oriz) + sin(orix)*sin(oriy)*sin(oriz)))/(cos(oriy)*cos(oriz)^2 + cos(oriy)*sin(oriz)^2)^2 - (cos(oriz)*(sin(oriy)*cos(oriz)^2 + sin(oriy)*sin(oriz)^2)*(cos(orix)*sin(oriz) - cos(oriz)*sin(orix)*sin(oriy)))/(cos(oriy)*cos(oriz)^2 + cos(oriy)*sin(oriz)^2)^2 + (cos(oriy)*cos(oriz)^2*sin(orix))/(cos(oriy)*cos(oriz)^2 + cos(oriy)*sin(oriz)^2) + (cos(oriy)*sin(orix)*sin(oriz)^2)/(cos(oriy)*cos(oriz)^2 + cos(oriy)*sin(oriz)^2)),0,0,0,0,0,0,0,0,0,0;0,0,0,- (sin(orix)*sin(oriz) + cos(orix)*cos(oriz)*sin(oriy))*(bay - amy + nay) - (cos(orix)*sin(oriz) - cos(oriz)*sin(orix)*sin(oriy))*(baz - amz + naz),cos(oriz)*sin(oriy)*(bax - amx + nax) - cos(orix)*cos(oriy)*cos(oriz)*(baz - amz + naz) - cos(oriy)*cos(oriz)*sin(orix)*(bay - amy + nay), (cos(orix)*cos(oriz) + sin(orix)*sin(oriy)*sin(oriz))*(bay - amy + nay) - (cos(oriz)*sin(orix) - cos(orix)*sin(oriy)*sin(oriz))*(baz - amz + naz) + cos(oriy)*sin(oriz)*(bax - amx + nax),0,0,0,0,0,0,0,0,0;0,0,0,(cos(oriz)*sin(orix) - cos(orix)*sin(oriy)*sin(oriz))*(bay - amy + nay) + (cos(orix)*cos(oriz) + sin(orix)*sin(oriy)*sin(oriz))*(baz - amz + naz),sin(oriy)*sin(oriz)*(bax - amx + nax) - cos(orix)*cos(oriy)*sin(oriz)*(baz - amz + naz) - cos(oriy)*sin(orix)*sin(oriz)*(bay - amy + nay), (cos(orix)*sin(oriz) - cos(oriz)*sin(orix)*sin(oriy))*(bay - amy + nay) - (sin(orix)*sin(oriz) + cos(orix)*cos(oriz)*sin(oriy))*(baz - amz + naz) - cos(oriy)*cos(oriz)*(bax - amx + nax),0,0,0,0,0,0,0,0,0;0,0,0,cos(oriy)*sin(orix)*(baz - amz + naz) - cos(orix)*cos(oriy)*(bay - amy + nay),cos(oriy)*(bax - amx + nax) + cos(orix)*sin(oriy)*(baz - amz + naz) + sin(orix)*sin(oriy)*(bay - amy + nay),0,0,0,0,0,0,0,0,0,0;0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;0,0,0,0,0,0,0,0,0,0,0,0,0,0,0];                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    
Ut=[0,0,0,0,0,0;0,0,0,0,0,0;0,0,0,0,0,0;sin(oriy) - (cos(oriy)*cos(oriz)^2*sin(oriy))/(cos(oriy)*cos(oriz)^2 + cos(oriy)*sin(oriz)^2) - (cos(oriy)*sin(oriy)*sin(oriz)^2)/(cos(oriy)*cos(oriz)^2 + cos(oriy)*sin(oriz)^2), (cos(oriz)*sin(oriy)*(cos(orix)*sin(oriz) - cos(oriz)*sin(orix)*sin(oriy)))/(cos(oriy)*cos(oriz)^2 + cos(oriy)*sin(oriz)^2) - cos(oriy)*sin(orix) - (sin(oriy)*sin(oriz)*(cos(orix)*cos(oriz) + sin(orix)*sin(oriy)*sin(oriz)))/(cos(oriy)*cos(oriz)^2 + cos(oriy)*sin(oriz)^2), (sin(oriy)*sin(oriz)*(cos(oriz)*sin(orix) - cos(orix)*sin(oriy)*sin(oriz)))/(cos(oriy)*cos(oriz)^2 + cos(oriy)*sin(oriz)^2) - (cos(oriz)*sin(oriy)*(sin(orix)*sin(oriz) + cos(orix)*cos(oriz)*sin(oriy)))/(cos(oriy)*cos(oriz)^2 + cos(oriy)*sin(oriz)^2) - cos(orix)*cos(oriy),0,0,0;0,- (cos(oriz)*(cos(orix)*cos(oriz) + sin(orix)*sin(oriy)*sin(oriz)))/(cos(oriz)^2 + sin(oriz)^2) - (sin(oriz)*(cos(orix)*sin(oriz) - cos(oriz)*sin(orix)*sin(oriy)))/(cos(oriz)^2 + sin(oriz)^2),(cos(oriz)*(cos(oriz)*sin(orix) - cos(orix)*sin(oriy)*sin(oriz)))/(cos(oriz)^2 + sin(oriz)^2) + (sin(oriz)*(sin(orix)*sin(oriz) + cos(orix)*cos(oriz)*sin(oriy)))/(cos(oriz)^2 + sin(oriz)^2),0,0,0;- (cos(oriy)*sin(oriz)^2)/(cos(oriy)*cos(oriz)^2 + cos(oriy)*sin(oriz)^2) - (cos(oriy)*cos(oriz)^2)/(cos(oriy)*cos(oriz)^2 + cos(oriy)*sin(oriz)^2),(cos(oriz)*(cos(orix)*sin(oriz) - cos(oriz)*sin(orix)*sin(oriy)))/(cos(oriy)*cos(oriz)^2 + cos(oriy)*sin(oriz)^2) - (sin(oriz)*(cos(orix)*cos(oriz) + sin(orix)*sin(oriy)*sin(oriz)))/(cos(oriy)*cos(oriz)^2 + cos(oriy)*sin(oriz)^2),(sin(oriz)*(cos(oriz)*sin(orix) - cos(orix)*sin(oriy)*sin(oriz)))/(cos(oriy)*cos(oriz)^2 + cos(oriy)*sin(oriz)^2) - (cos(oriz)*(sin(orix)*sin(oriz) + cos(orix)*cos(oriz)*sin(oriy)))/(cos(oriy)*cos(oriz)^2 + cos(oriy)*sin(oriz)^2),0,0,0;0,0,0, -cos(oriy)*cos(oriz),cos(orix)*sin(oriz) - cos(oriz)*sin(orix)*sin(oriy),- sin(orix)*sin(oriz) - cos(orix)*cos(oriz)*sin(oriy);0,0,0,-cos(oriy)*sin(oriz), - cos(orix)*cos(oriz) - sin(orix)*sin(oriy)*sin(oriz),cos(oriz)*sin(orix) - cos(orix)*sin(oriy)*sin(oriz);0,0,0,sin(oriy),-cos(oriy)*sin(orix),-cos(orix)*cos(oriy);1,0,0,0,0,0;0,1,0,0,0,0;0,0,1,0,0,0;0,0,0,1,0,0;0,0,0,0,1,0;0,0,0,0,0,1];
Ft=eye(15)+dt*At;
Vt=Ut;
Q=diag([nbgx nbgy nbgz nbax nbay nbaz]);
Qd=Q*dt;
f=x_dot;
uEst=uPrev+dt*f;
covarEst=Ft*covarPrev*transpose(Ft)+Vt*Qd*transpose(Vt);
end


