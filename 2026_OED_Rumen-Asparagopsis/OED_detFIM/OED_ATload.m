
% Muñoz-Tamayo 2026

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Loading the data 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


global tDMI yDMI 

Ny= 7; 

% Initial conditions of the state variables 

Cinit = [  1.9552  0.93422  0.22752  0.0011434  0.0001231  0.028303  0.0079153  0.012547  0.013099  0.00064129  1.3357e-06  0.0002533  0.0043116  0.0031899  0.00030976  0.0015073  4.9056e-06  0.00063912  0  ];


% Setting of the DMI 

DMtotal = 11.25; %total dry matter ingested in one day, g
k = 0.015; %the intake rate constant, 1/h

%ni = 1; %number of feed distributions in one day
%ww = [1]; %Distribution of DMtotal between meals

%ni = 2; %number of feed distributions in one day
%ww = [0.50,0.50]; %Distribution of DMtotal between meals

ni = 2; %number of feed distributions in one day
ww = [0.70,0.30]; %Distribution of DMtotal between meals

%ni = 4; %number of feed distributions in one day
%ww = [0.25,0.25,0.25,0.25]; %Distribution of DMtotal between meals


DMIdyn = dynamicIntake(DMtotal,ww,ni,k); 
tDMI = DMIdyn(:,1); %time of measurements of DMI, h 
yDMI = DMIdyn(:,2); %measurements of DMI, g/h 
 
% % plotting the DMI 
% figure
% plot(tDMI,yDMI,'linewidth', 2);
% set(gca,'Fontsize',18);
% 	 xlabel(' Time (h) ', 'fontsize',16);
% 	 ylabel('DMI (g/h)');
