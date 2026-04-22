% Muñoz-Tamayo 2026
% program to perform the optimization, maximization of the determinant of
% the FIM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Optimization 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all


alt_optim= 1; 


load detFIMparamOLD.mat
load detFIMparam.mat

% Estimates with IS = 1/(5e-4)^2

load detFIMparam_optim_n1_5.mat     % 1 feed distribution, det(F) = 3.5087e+20
load detFIMparam_optim_n2a_5.mat    % 2 feed distributions at 50%, det(F) = 2.7881e+18
load detFIMparam_optim_n2b_5.mat    % 2 feed distributions at 70%, 30% det(F) = 1.3714e+19
load detFIMparam_optim_n4_5.mat     % 4 feed distributions at 25%, det(F) = 1.1883e+19

detFIMparam0 = [(24/5)*ones(1,5)-0.5  0.0034]; % equidistant sampling

% Setting the inital guess of the estimates 
     
     paramMi = detFIMparam;  
     %paramMi = detFIMparamOLD 
     %paramMi = detFIMparam0;
     %paramMi = detFIMparam_optim_n1_5;
     %paramMi = detFIMparam_optim_n2a_5;
     %paramMi = detFIMparam_optim_n2b_5;
     %paramMi = detFIMparam_optim_n4_5;
     
     
     paramMi = log(paramMi);
     

tic;
disp('                                       ');
	 % disp('   Type 1 i you are using normal scale, or 2 if it is in log scale : '); 
	 % olog = input('    '); 
	 olog = 1; 

clc ; 
	 if olog==1
	 LBp = paramMi*0.9;
	 UBp = paramMi*1.1;
	 else 
	 LBp = log(exp(paramMi*0.8));
	 UBp = log(exp(paramMi*1.2));
	 end 

	 N = 1; % Number of runs of the optimization algorithm 
	 options1 = optimset('LargeScale','off', 'Display','iter', 'Diagnostics','on', 'MaxIter', 100);
	 options2 = optimset('LargeScale','on', 'GradObj','on','Hessian', 'on', 'Display','iter', 'Diagnostics','on');

	 if alt_optim ==1 
	 options = options1; 
	 else 
	 options = options2; 
	 end 
	 clear i; 
disp('                                       ');
disp('           *********************************** ');
disp('                  Finding the estimates');
disp('           *********************************** ');
	 
Np = length(paramMi);

for i = 1:N 
		  aleat = rand(1,Np);
 		  if i==1
 		  paramo(1,:) = paramMi; 
 		  else 
 		   for j=1:Np  
 		    paramo(i,j) = LBp(j)+aleat(j)*(UBp(j)-LBp(j));  
 		  end
       
          
 		  end 
 	 %[x(i,:),fval(i), exitflag(i), output,grad,hessian] = fminunc(@det_FIMcost,paramo(i,:),options);
	 [x(i,:),fval(i), exitflag(i), output] = fminsearch(@det_FIMcost,paramo(i,:),options);
disp('                                       ');
	 fvalues(i) = fval(i);
	 estimates(i,:) = x(i,:);
		  end 
 
	 F = sort(fvalues);
	 minfval = min(fvalues); 
	 indexfval = find(fvalues==minfval); 
	 teta = estimates(indexfval,:);
     %teta(5) =  52.0946; % setting the initial time 
	 Jc = minfval; 
	 teta_estimates = teta;
     teta_estimates = exp(teta); % For the case of change in the parameterization 
	 zero_teta = length(find(teta_estimates>=0)); 
	 J_cost = Jc; 

tiempo = toc;
detFIMparamOLD = paramMi; 
detFIMparam = teta_estimates; 

save detFIMparam.mat detFIMparam
save detFIMparamOLD.mat detFIMparamOLD

deltaT = 0.5; % setting an interval of at least 30 min between sampling times

ts(1) = 0; 
clear k 
ns = length(detFIMparam); 

for k=2:ns
ts(k) = ts(k-1) + detFIMparam(k-1) + deltaT;
end
ts(end) = round(ts(end));
ts
