% Muñoz-Tamayo 2026
% Dynamic model of the dry matter intake (g/h) for one day

function [x] = dynamicIntake(DMtotal,ww,ni,k)

%DMtotal: total dry matter ingested in one day, g
%ww: distribution of DMtotal between meals, %
%ni: number of feed distributions in one day
%k: intake rate constant, 1/h

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Defining the dry matter intake (g) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

time = 0:1/60:24; % time in hour
tm = mod(time,(time(end)/ni)); % applies modulo function to work with nm symetric responses

DM=[];
for i = 1:ni
    DM_max=ww(i)*DMtotal;
    indexT = find(tm==0); 
    DM=[DM,DM_max*exp(-k*60*tm(indexT(i):indexT(i+1)-1))];
end
DM((time(end)*60)+1)=ww(1)*DMtotal;

%Computation of DMI in g/h (DMI = d/dt DM)
DMI(1) = 0; %initial condition of DMI
for j = 2:length(DM)
    DMI(j) = (DM(j-1)-DM(j))/(time(2)-time(1)); 
    DMI(j) = max(0,DMI(j)); 
end

x = [time' DMI'];

end 
