function [uCurr,covar_curr] = upd_step(z_t,covarEst,uEst)
%z_t is the measurement
%covarEst and uEst are the predicted covariance and mean respectively
%uCurr and covar_curr are the updated mean and covariance respectively
%%z=[1 0 0 0 0 0 0 0 0 0 0 0 0 0 0;0 1 0 0 0 0 0 0 0 0 0 0 0 0 0;0 0 1 0 0 0 0 0 0 0 0 0 0 0 0;0 0 0 1 0 0 0 0 0 0 0 0 0 0 0;0 0 0 0 1 0 0 0 0 0 0 0 0 0 0;0 0 0 0 0 1 0 0 0 0 0 0 0 0 0]*uPrev+[npx;npy;npz;nroll;npitch;nyaw];
ct=[1 0 0 0 0 0 0 0 0 0 0 0 0 0 0;0 1 0 0 0 0 0 0 0 0 0 0 0 0 0;0 0 1 0 0 0 0 0 0 0 0 0 0 0 0;0 0 0 1 0 0 0 0 0 0 0 0 0 0 0;0 0 0 0 1 0 0 0 0 0 0 0 0 0 0;0 0 0 0 0 1 0 0 0 0 0 0 0 0 0];
Kt=covarEst*transpose(ct)*inv(ct*covarEst*transpose(ct)+eye(6)*(0.0001*eye(6))*transpose(eye(6)));
uCurr=uEst+Kt*(z_t-ct*uEst);
covar_curr=covarEst-Kt*ct*covarEst;
end