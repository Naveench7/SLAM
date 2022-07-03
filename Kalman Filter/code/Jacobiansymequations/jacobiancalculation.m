%jacobian calculation using symbolic math and results are used in pred_step
syms orix oriy oriz vx vy vz bgx bgy bgz bax bay baz x y z
syms G omegamx omegamy omegamz ngx ngy ngz g gx gy gz R amx amy amz nax nay naz nbgx nbgy nbgz nbax nbay nbaz
V=[vx;vy;vz];
g=[gx;gy;gz];
G=[0 -sin(oriz) cos(oriz)*cos(oriy);0 cos(oriz) sin(oriz)*cos(oriy);1 0 -sin(oriy)];
Gi=inv(G);
simplify(Gi);
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
x_dot=[V;Gi*R*(Omegam-Bg-Ng);g+(R*(Am-Ba-Na));Ng;Na];
xs=[x;y;z;orix;oriy;oriz;vx;vy;vz;nbgx;nbgy;nbgz;nbax;nbay;nbaz];
at=jacobian(x_dot,xs)
n=[ngx;ngy;ngz;nax;nay;naz];
ut=jacobian(x_dot,n)

