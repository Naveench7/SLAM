 function [covarEst,uEst] = pred_step(uPrev,covarPrev,angVel,acc,dt)
%% BEFORE RUNNING THE CODE CHANGE NAME TO pred_step
    %% Parameter Definition
    % uPrev - is the mean of the prev state
    %covarPrev - covar of the prev state
    %angVel - angular velocity input at the time step
    %acc - acceleration at the timestep
    %dt - difference in time 
     g= [0; 0; -9.81];
    wm = [angVel(1); angVel(2); angVel(3)];
    am = [acc(1); acc(2); acc(3)];
    alpha=0.001;
    kappa=1;
    Qt=1;
    beta=2;
    
    %computing sigma points
    n_dash=27;
    lambda = ((alpha)^2)*(n_dash+kappa)-n_dash;
    u_aug = [uPrev;zeros(12,1)];
    covar_aug = [covarPrev,zeros(15,12);zeros(12,15),eye(12)*Qt];
    X_p = [u_aug];
    ch_covar = chol(covar_aug,"lower");
    X_p = (X_p + [zeros(n_dash,1) sqrt(n_dash+lambda)*ch_covar -sqrt(n_dash+lambda)*ch_covar]);
   
   
    %propagating sigma points through the non linear function
    X_t=[];
    for j=1:2*n_dash+1
        x_aug=X_p(1:15,j);
        n_aug=X_p(16:27,j);
        v=[x_aug(7);x_aug(8);x_aug(9)];
        bg=[x_aug(10);x_aug(11);x_aug(12)];
        ba=[x_aug(13);x_aug(14);x_aug(15)];
        ng=[n_aug(1);n_aug(2);n_aug(3)];
        na=[n_aug(4);n_aug(5);n_aug(6)];
        nbg=[n_aug(7);n_aug(8);n_aug(9)];
        nba=[n_aug(10);n_aug(11);n_aug(12)];
        roll=x_aug(4);
        pitch=x_aug(5);
        yaw=x_aug(6);
        G =[0,-sin(yaw),cos(pitch)*cos(yaw);
            0,cos(yaw),cos(pitch)*sin(yaw);
            1,0,-sin(pitch)];
        R =[cos(pitch)*cos(yaw), cos(yaw)*sin(pitch)*sin(roll) - cos(roll)*sin(yaw), sin(roll)*sin(yaw) + cos(roll)*cos(yaw)*sin(pitch); 
            cos(pitch)*sin(yaw), cos(roll)*cos(yaw) + sin(pitch)*sin(roll)*sin(yaw), cos(roll)*sin(pitch)*sin(yaw) - cos(yaw)*sin(roll); 
            -sin(pitch), cos(pitch)*sin(roll), cos(pitch)*cos(roll)];
        Gi=inv(G);
        x_t=x_aug+[v;Gi*R*(wm-bg-ng);g+(R*(am-ba-na));nbg;nba]*dt;
        X_t=[X_t,x_t];
    end
    
   %computing the predicted mean,predicted covariance of the
    %meaasurement,and predicted cross-covariance 
    
    W_0m = lambda/(n_dash+lambda);
    W_im = 1/(2*(n_dash+lambda));
    for m= 1:2*n_dash+1
       if m==1
           initial=W_0m*X_t(:,m);
           pinitial= initial;
       else
           remaining=W_im*X_t(:,m);
           pinitial= pinitial+remaining;
       end
    
    end
    uEst=pinitial;

    W_0c= (lambda/(n_dash+lambda))+(1-(alpha^2)+ beta);
    W_ic=1/(2*(n_dash+lambda));
    for n= 1:2*n_dash+1
       if n==1
            initial_covar=W_0c*(X_t(:,n)-uEst)*(X_t(:,n)-uEst)';
            p_cinitial=initial_covar;
       else
           remaining_covar=W_ic*(X_t(:,n)-uEst)*(X_t(:,n)-uEst)';
            p_cinitial=p_cinitial+remaining_covar;
       end
    end
   covarEst=p_cinitial;
   
end

