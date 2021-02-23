%% Seasonal to Mean Annual Correction Method 
%% About
%
% *SAT*: This MATLAB function converts seasonal to mean annual
% temperature. The seasonality of the record must be prescribed. If 
% seasonality is unknown use SAT_v2.m
%For more information read reference below.
%
% The matlab script daily_insolation.m to compute daily average insolation 
% written by Huybers and Eisenman (2006) is required to run this function 
% (http://eisenman.ucsd.edu/code.html).
%
% Written by Samantha Bova (samantha.bova@rutgers.edu), July, 2020. Updated
% February 2021.
%
% CITATION: Bova, S., Rosenthal, Y., Liu, Z. et al. Seasonal origin of 
% the thermal maxima at the Holocene and the last interglacial. Nature 589, 
% 548?553 (2021). https://doi.org/10.1038/s41586-020-03155-x
%% Input
% SST data should be input as a two column matrix with the following format:
% * Column 1 - time (kyrs BP only)
% * Column 2 - SSTsn (°C)
% Remove all NaNs as they will cause an error.
% To run type the following in MATLAB command line:  
% >  [core_ma]=SAT_v1(core,latitude,SSTsens,seas_start,seas_end);
%
% replace core with name of two column matrix with your SSTsn data 
% replace latitude with latitude for the core site in decimal degrees
% replace SSTsens with SST sensitivity
% replace seas_start with first day of year of season represented by
% record
% replace seas_end with last day of year of season represented by record
%% Output
% core_ma = [age (ky), MASST (°C), error (1 standard error)];
%%
function [core_ma]=SAT_seas_prec_v1(core,latitude,SSTsens,seas_start,seas_end)

latmin=latitude-0.5;
latmax=latitude+0.5;

core_ma2=[];

% MA Insolation forcing
[kyear, lat, day]=meshgrid(0:1:150,latmin:0.1:latmax,1:365);
     MA_Inso=daily_insolation(kyear,lat,day);
     MA_Inso=mean(MA_Inso,1);
     MA_Inso=mean(MA_Inso,3);
     MA_Inso=squeeze(MA_Inso);
     MAinsoi=interp1(0:1:150,MA_Inso,core(:,1));
     MAanomi=MAinsoi-MAinsoi(1);

% Records identified as seasonal will are adjusted to mean annual

Seasonality_start_end=[seas_start,seas_end]

[kyear, lat, day]=meshgrid(0:1:150,latmin:0.1:latmax,Seasonality_start_end(1):Seasonality_start_end(2));
Seas_Inso=daily_insolation(kyear,lat,day);
Seas_Inso=mean(Seas_Inso,1);
Seas_Inso=mean(Seas_Inso,3);
Seas_Inso=squeeze(Seas_Inso);
seasinsoi=interp1(0:1:150,Seas_Inso,core(:,1));
seasanomi=seasinsoi-seasinsoi(1);
delta_I=seasanomi-MAanomi;

for i=1:length(SSTsens)
% Apply linear SST Sensitvity to seasonal insolation
SSTSensitivity=SSTsens(i);

% Calculate SST adjustment to convert seasonal to mean annual SST
MA_anom=(MA_Inso(1:151)-MA_Inso(1));
Seas_anom=(Seas_Inso(1:151)-Seas_Inso(1));
SSTadjustment=(Seas_anom'-MA_anom').*(SSTSensitivity);
SSTadjustmenti=interp1(0:150,SSTadjustment,core(:,1),'linear','extrap');

core_MA=core(:,2)-SSTadjustmenti;
core_ma2=[core_ma2,core_MA];
end

core_ma=[core(:,1),mean(core_ma2,2),(std(core_ma2,0,2)./sqrt(2))];

figure(3)
plot(core_ma(:,1),core_ma(:,2),'-o','LineWidth',3)
hold on
plot(core(:,1),core(:,2),'-o','LineWidth',3)
set(gca,'FontSize',12)
ylabel('SST Anomaly (°C)')
set(gca,'xdir','reverse') 
xlim([0 12])
xlabel('Age(ka)')
title('Holocene')
set(gca,'xdir','reverse') 
legend('Mean Annual','Original unadjusted data')
end
