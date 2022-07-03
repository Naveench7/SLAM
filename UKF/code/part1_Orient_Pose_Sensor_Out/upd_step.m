function [uCurr,covar_curr] = upd_step(z_t,covarEst,uEst)
%% BEFORE RUNNING THE CODE CHANGE NAME TO upd_step
    %% Parameter Definition
    %z_t - is the sensor data at the time step
    %covarEst - estimated covar of the  state
    %uEst - estimated mean of the state
    c=eye(6,15);
    Kt=covarEst*c'*inv(c*covarEst*c'+0.01*eye(6));
    covar_curr=covarEst-Kt*c*covarEst;
    uCurr=uEst+Kt*(z_t-c*uEst);  
end

