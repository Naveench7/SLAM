function [uCurr,covar_curr] = upd_step(z_t,covarEst,uEst)
%% BEFORE RUNNING THE CODE CHANGE NAME TO upd_step
    %% Parameter Definition
    %z_t - is the sensor data at the time step
    %covarEst - estimated covar of the  state
    %uEst - estimated mean of the state
    
    %initialization
    alpha = 0.001; 
    kappa = 1; 
    beta = 2; 
    
    %%Using rotation and translantion of camera w.r.t body from proj 2 to
   %%convert into body w.r.t camera
    t_cb = [-0.04;0;-0.03];
    R_cb = [1.414/2,-1.414/2,0;-1.414/2,-1.414/2,0;0,0,-1]; 
    R_bc = (R_cb)';
    H_cb = [R_cb,t_cb;0,0,0,1];
    Hb_c = inv(H_cb);
    
    %computing Sigma Points
    n_ddash = 15; 
    lambda_ddash =(alpha^2 * (n_ddash + kappa)) - n_ddash;
    u_auge = uEst;
    covar_auge = covarEst;
    cho_covar = chol(covar_auge,"lower");
    Q =  0.04 * eye(3); 
    X_i=[u_auge];
    X_i=(X_i + [zeros(n_ddash,1) sqrt(n_ddash+lambda_ddash)*cho_covar -sqrt(n_ddash+lambda_ddash)*cho_covar]);
   
    
    %propagating sigma points through the non linear function
    Z_t=[];
    for zu = 1 :(2*n_ddash+1)
       x_auge=X_i(1:15,zu);
       roll = x_auge(4);
       pitch = x_auge(5);
       yaw = x_auge(6);
       R=[cos(pitch)*cos(yaw), cos(yaw)*sin(pitch)*sin(roll) - cos(roll)*sin(yaw), sin(roll)*sin(yaw) + cos(roll)*cos(yaw)*sin(pitch);
           cos(pitch)*sin(yaw), cos(roll)*cos(yaw) + sin(pitch)*sin(roll)*sin(yaw), cos(roll)*sin(pitch)*sin(yaw) - cos(yaw)*sin(roll);
           -sin(pitch), cos(pitch)*sin(roll), cos(pitch)*cos(roll)];
       zt = (R_cb*R'*x_auge(7:9))-(R_cb*skew(Hb_c(1:3,4))*R_bc*z_t(4:6));
       Z_t=[Z_t zt];   
   end

   %computing the predicted mean,predicted covariance of the
    %meaasurement,and predicted cross-covariance
  
  W0_mz = lambda_ddash/(n_ddash + lambda_ddash);
  Wi_mz = 1/(2 * (n_ddash + lambda_ddash));

  for zu=1:(2*n_ddash + 1)
      if zu==1
          initial=W0_mz*Z_t(:,zu);
          z =initial;
      else
          remaining=Wi_mz * Z_t(:,zu);
          z = z +remaining;
      end
  end
  
  
  W0_cz = (lambda_ddash/(n_ddash + lambda_ddash)) + (1 - alpha^2 - beta);
  Wi_cz = 1/(2 * (n_ddash + lambda_ddash));

  for zu=1:(2*n_ddash + 1)
      if zu == 1
          cinitial=W0_cz*(X_i(1:15 ,zu)-uEst)*(Z_t(:, zu)-z)';
          Ct = cinitial;
          sinitial=W0_cz*(Z_t(:,zu)-z)*(Z_t(:,zu)-z)';
          St = sinitial;
      else
          cremaining=(Wi_cz*(X_i(1:15 ,zu)-uEst)*(Z_t(:,zu)-z)');
          Ct = Ct + cremaining;
          sremaining=(Wi_cz*(Z_t(:, zu)-z)*(Z_t(:,zu)-z)');
          St = St + sremaining;

      end

  end

  St = St+Q;

  %Computing the filter gain and the filtered state mean and
  % covariance,conditional to the measurement
  Kt = Ct/St;
  uCurr = uEst+(Kt*(z_t(1:3)-z));
  covar_curr =covarEst-(Kt*St*transpose(Kt));

end

    

