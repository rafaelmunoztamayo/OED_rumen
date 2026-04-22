% Muñoz-Tamayo 2026
% program to calculate de determinant of the FIM 

clear all

global wBrOptim
OED_ATload; 

load detFIMparam.mat
detFIMparam0 = [(24/5)*ones(1,5)-0.5  0.0034]; % equidistant sampling

% Estimates with IS = 1/(5e-4)^2

load detFIMparam_optim_n1_5.mat     % 1 feed distribution, det(F) = 3.5087e+20
load detFIMparam_optim_n2a_5.mat    % 2 feed distributions at 50%, det(F) = 2.7881e+18
load detFIMparam_optim_n2b_5.mat    % 2 feed distributions at 70%, 30% det(F) = 1.3714e+19
load detFIMparam_optim_n4_5.mat     % 4 feed distributions at 25%, det(F) = 1.1883e+19

% standard deviations with equidistant sampling from detFIMparam0 
load sDeqn1.mat         % 1 feed distribution, det(F) = 4.0332e+19
load sDeqn2a.mat        % 2 feed distributions at 50%, det(F) = 5.6142e+17
load sDeqn2b.mat        % 2 feed distributions at 70%, 30% det(F) = 1.0377e+18
load sDeqn4.mat         % 4 feed distributions at 25%, det(F) = 2.7790e+18

load sDn1_5.mat       % 1 feed distribution, det(F) = 3.5087e+20
load sDn2a_5.mat      % 2 feed distributions at 50%, det(F) = 2.7881e+18
load sDn2b_5.mat      % 2 feed distributions at 70%, 30% det(F) = 1.3714e+19
load sDn4_5.mat       % 4 feed distributions at 25%, det(F) = 1.1883e+19

% determinants of the FIM
dFeqn1  = 4.0332e+19;
dFeqn2a = 5.6142e+17; 
dFeqn2b = 1.0377e+18; 
dFeqn4  = 2.7790e+18; 

dFoptimn1  = 3.5087e+20; 
dFoptimn2a = 2.7881e+18; 
dFoptimn2b = 1.3714e+19; 
dFoptimn4  = 1.1883e+19; 

dFeqn  = [dFeqn1 dFeqn2a dFeqn2b dFeqn4 ];
dF_opt = [dFoptimn1 dFoptimn2a dFoptimn2b dFoptimn4 ]; 

% ratio of det(F)
clear j 
for j=1:length(dF_opt)
ratio_detF(j) = [dF_opt(j)./dFeqn(j)];
end

% ratio of sD

ratio_sD(1,:) = [100*sDn1_5./sDeqn1];  
ratio_sD(2,:) = [100*sDn2a_5./sDeqn2a];  
ratio_sD(3,:) = [100*sDn2b_5./sDeqn2b];  
ratio_sD(4,:) = [100*sDn4_5./sDeqn4];  

ratio_sD_best = [100*sDn1_5./sDeqn2b];

Np= 5;  % number of parameters
Nx= 19; % number of state variables

% selecting the parameters 

     %paramt = detFIMparam; 
     %paramt = detFIMparam0;
          
     %paramt = detFIMparam_optim_n1_5; 
     %paramt = detFIMparam_optim_n2a_5; 
     paramt = detFIMparam_optim_n2b_5; 
     %paramt = detFIMparam_optim_n4_5;
     
     wBrOptim = paramt(end); % including the dose of bromoform in the OED



% calculation of the sampling times 
deltaT = 0.5; % setting an interval of at least 30 min between sampling times
ts(1) = 0; 
clear k 

%ns = 10; % number of samples
%ns = 5; % number of samples
ns = length(paramt) ; 

for k=2:ns
ts(k) = ts(k-1) + paramt(k-1) + deltaT;
end

ts(end) = round(ts(end));

ts

final_deltas = [detFIMparam_optim_n1_5; detFIMparam_optim_n2a_5; detFIMparam_optim_n2b_5; detFIMparam_optim_n4_5]; 
tTs = zeros(4,1);
wBrfinal = [detFIMparam_optim_n1_5(end); detFIMparam_optim_n2a_5(end); detFIMparam_optim_n2b_5(end); detFIMparam_optim_n4_5(end)]; 
for j=1:4
for k=2:ns
tTs(j,k) = tTs(j,k-1) + final_deltas(j,k-1) + deltaT;
tTs(j,k) = round(tTs(j,k)*10)/10;
tTs(j,k) = min(tTs(j,k),24);
end
end


x_br_AT   = 6.84/1000;          % Content of bromoform in Asparagopsis taxiformis, g/g
wATfinal = wBrfinal./x_br_AT; % Fraction of AT

final_params = [wATfinal tTs(:,2:end)];


ceros=zeros(1,Np*Nx); 

Cin = [Cinit  ceros]; 
[ti,Ci]=ode15s('OED_ATse',[ts],[Cin]);

Xm=Ci(2:end,1:Nx); 
%Xmp=Ci(1:end,1:Nx); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Computation of outputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear k 
for k=1:length(ti)-1
	 Ym(:,k) = OED_ATout(Xm(k,:));
end 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Computation of sensitivities
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sensitivities of the state
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	 s_1=(Ci(2:end,20:38))';
	 s_2=(Ci(2:end,39:57))';
	 s_3=(Ci(2:end,58:76))';
	 s_4=(Ci(2:end,77:95))';
	 s_5=(Ci(2:end,96:114))';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sensitivities of the output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 clear k
 
for k=1:length(ti)-1
[dymds_1(:,k), dymds_2(:,k), dymds_3(:,k), dymds_4(:,k), dymds_5(:,k),] = OED_ATsy(Xm(k,:),s_1(:,k),s_2(:,k),s_3(:,k),s_4(:,k),s_5(:,k));
end

[rowy,coly]=size(Ym);
nt=length(ti);

    Sigma = 5e-4^2; %7.7e-4^2; % estimated as 10% of the average mean of the measurements 5e-4 is between of 7.7 and mean(sdY)
    IS = inv(Sigma); 
      
    % Calculation of the Fisher Information Matrix FIM 
    
    F = zeros(Np,Np);
	F = zeros(5,5);
	  for  i=1:coly 
	    dymdp = [	 dymds_1(:,i) 	 dymds_2(:,i) 	 dymds_3(:,i) 	 dymds_4(:,i) 	 dymds_5(:,i) ];
	    F  = dymdp'*IS*dymdp +F;
	    
      end 
 
      dF = -det(F); % determinant of the FIM 
     
      if ts(end)>24
          dF = 0; 
      end 
      
      if (ts(end)-ts(end-1))<deltaT
          dF = 0; 
      end 

 dF     
% G = -G;
 P = inv(F); % Covariance matrix of the estimator 
 sP = P.^0.5;
 sD = diag(sP); %standard deviation of the estimates 


% plotting the model reponse 

% figure; 
% 
% [ti2,Ci2]=ode15s('OED_ATse',[ts(1) ts(end)],[Cin]);
% clear k 
% Xmr=Ci2(1:end,1:Nx); 
% for k=1:length(ti2)
% 	 Ymm(:,k) = OED_ATout(Xmr(k,:));
% end 
% 
% clear k 
% Xmr=Ci(1:end,1:Nx); 
% for k=1:length(ti)
% 	 Ymr(:,k) = OED_ATout(Xmr(k,:));
% end
% subplot(241), plot(ti,Ymr(1,:),'o','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',8);
% hold on 
% plot(ti2,Ymm(1,:),'linewidth', 2);
% set(gca,'fontsize',20);
% 	 xlabel(' Time (h) ', 'fontsize',20);
% 	 ylabel(' [acetate] mol/L ', 'fontsize',20);
% subplot(242), plot(ti,Ymr(2,:),'o','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',8);
% hold on 
% plot(ti2,Ymm(2,:),'linewidth', 2);
% set(gca,'fontsize',20);
% 	 xlabel(' Time (h) ', 'fontsize',20);
% 	 ylabel(' [butyrate] mol/L ', 'fontsize',20);
% subplot(243), plot(ti,Ymr(3,:),'o','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',8);
% hold on 
% plot(ti2,Ymm(3,:),'linewidth', 2);
% set(gca,'fontsize',20);
% 	 xlabel(' Time (h) ', 'fontsize',20);
% 	 ylabel(' [propionate] mol/L ', 'fontsize',20);
%      
% subplot(244),plot(ti,Ymr(7,:),'o','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',8);
% hold on 
% plot(ti2,Ymm(7,:),'linewidth', 2);
% set(gca,'fontsize',20);
% 	 xlabel(' Time (h) ','fontsize',20);
% 	 ylabel(' [bromoform] g/L ','fontsize',20);
% 
%      subplot(245), plot(ti,Ymr(4,:),'o','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',8);
% hold on 
% plot(ti2,Ymm(4,:),'linewidth', 2);
% set(gca,'fontsize',20);
% 	 xlabel(' Time (h) ', 'fontsize',20);
% 	 ylabel(' qCO_2 (mol/h) ', 'fontsize',20);
% subplot(246), plot(ti,Ymr(5,:),'o','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',8);
% hold on 
% plot(ti2,Ymm(5,:),'linewidth', 2);
% set(gca,'fontsize',20);
% 	 xlabel(' Time (h) ', 'fontsize',20);
% 	 ylabel(' qH_2 (mol/h)');
% subplot(247), plot(ti,Ymr(6,:),'o','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',8);
% hold on 
% plot(ti2,Ymm(6,:),'linewidth', 2);
% set(gca,'fontsize',20);
% 	 xlabel(' Time (h) ', 'fontsize',20);
% 	 ylabel(' qCH_4 (mol/h)');
% 
