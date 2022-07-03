function [uCurr,covar_curr] = upd_step(z_t,covarEst,uEst)
%z_t is the measurement
%covarEst and uEst are the predicted covariance and mean respectively
%uCurr and covar_curr are the updated mean and covariance respectively
%ct=[0 0 0 0 0 0 1 0 0 0 0 0 0 0 0;0 0 0 0 0 0 0 1 0 0 0 0 0 0 0;0 0 0 0 0 0 0 0 1 0 0 0 0 0 0];
ct=[0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;0 0 0 0 0 0 1 0 0 0 0 0 0 0 0;0 0 0 0 0 0 0 1 0 0 0 0 0 0 0;0 0 0 0 0 0 0 0 1 0 0 0 0 0 0];
Kt=covarEst*transpose(ct)*inv(ct*covarEst*transpose(ct)+eye(9)*(0.00001*eye(9))*transpose(eye(9)));
uCurr=uEst+Kt*([0;0;0;0;0;0;z_t]-ct*uEst);
covar_curr=covarEst-Kt*ct*covarEst;
end