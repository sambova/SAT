%% Seasonal to Mean Annual Correction Method 
%% About
%
% *SAT*: This MATLAB function assesses the seasonality of temperature-sensitive
% proxy records from LIG records and converts seasonal to mean annual
% temperature. For more information read reference below.
%
% The matlab script daily_insolation.m to compute daily average insolation writtent by Huybers
% and Eisenman (2006) is required to run this function (http://eisenman.ucsd.edu/code.html).
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
% >  [core_ma,SST_sensitivity,Start_day,End_day,Linear_r2]=SAT_v1(core*,latitude*);
%
% * replace core with name of two column matrix with your SSTsn data and
% latitude with latitude for the core site in decimal degrees
%% Output
% core_ma = [age (ky), MASST (°C), error (1 standard
% error)];
% SSTsens = [SST sensitivity, 1 std. dev.];
% Start_day = day of year at start of 30 day window of best fit
% End_day = day of year at end of 30 day window of best fit
% Linear R2 = average r2 of linear regression between LIG SSTs and 30 day window of best
% fit
%%
function [core_ma,SST_sensitivity,Start_day,End_day,Linear_r2]=SAT_v1_check(core,latitude)

latmin=latitude-0.5;
latmax=latitude+0.5;

% Find SST data from LIG (115-127 ka)
indlig=find(core(:,1)>115 & core(:,1)<127);
ind=length(indlig);
ssti=core(indlig,2);

% MA Insolation forcing
[kyear, lat, day]=meshgrid(0:1:150,latmin:0.1:latmax,1:365);
     MA_Inso=daily_insolation(kyear,lat,day);
     MA_Inso=mean(MA_Inso,1);
     MA_Inso=mean(MA_Inso,3);
     MA_Inso=squeeze(MA_Inso);
     MAinsoi=interp1(0:1:150,MA_Inso,core(indlig,1));
     MAanomi=MAinsoi-MAinsoi(1);

core_ma2=[];
SSTsens=[];
start_day=[];
end_day=[];

for ii=1:ind;
    indlig=find(core(:,1)>115 & core(:,1)<127);
    indlig(ii)=[];
    ssti=core(indlig,2);
% Test for month of best fit using LIG SSTs
r2s=[];
slopes=[];
%days=1:730;
w=30;
i=0:730;
for g=1:10:400;
    days=i(g):i(g+30);
    [kyear, lat, day]=meshgrid(0:1:150,latmin:0.1:latmax,days);
     Seas_Inso=daily_insolation(kyear,lat,day);
     Seas_Inso=mean(Seas_Inso,1);
     Seas_Inso=mean(Seas_Inso,3);
     Seas_Inso=squeeze(Seas_Inso);
     seasinsoi=interp1(0:1:150,Seas_Inso,core(indlig,1));
     MAinsoi=interp1(0:1:150,MA_Inso,core(indlig,1));
     MAanomi=MAinsoi-MAinsoi(1);
     seasanomi=seasinsoi-seasinsoi(1);
     delta_I=seasanomi-MAanomi;
     p=polyfit(delta_I,ssti,1);
     slope=p(1);
     slopes=[slopes,slope];
     mdl=fitlm(ssti,delta_I,'linear');
     r2=mdl.Rsquared.Ordinary;
     r2s=[r2s,r2];
end
%%
ind=find(slopes>=0);
Seasmax_r2=max(r2s(ind));
ind2=find(r2s==Seasmax_r2);
seas_start=(ind2-1)*10;
seas_end=seas_start+30;

start_day=[start_day,seas_start];
end_day=[end_day,seas_end];

% if record shows best fit with mean annual temperature (12 months) no
% adjustments are applied and  original data is output
     map=polyfit(MAanomi,ssti,1);
     maslope=map(1);
     mdl=fitlm(ssti,MAanomi,'linear');
     MAr2=mdl.Rsquared.Ordinary;
     
if maslope > 0 && MAr2 >= Seasmax_r2 
    core_MA=core(:,2);
    SST_Sensitivity=1./maslope;
    MAr2=mdl.Rsquared.Ordinary;
    start_day=[];
    end_day=[];
    disp('Record reflects Mean Annual Temperature')
    
% Records identified as seasonal will are adjusted to mean annual
elseif maslope < 0 || MAr2 < Seasmax_r2  
Seasonality_start_end=[seas_start,seas_end];

[kyear, lat, day]=meshgrid(0:1:150,latmin:0.1:latmax,Seasonality_start_end(1):Seasonality_start_end(2));
Seas_Inso=daily_insolation(kyear,lat,day);
Seas_Inso=mean(Seas_Inso,1);
Seas_Inso=mean(Seas_Inso,3);
Seas_Inso=squeeze(Seas_Inso);
seasinsoi=interp1(0:1:150,Seas_Inso,core(indlig,1));
seasanomi=seasinsoi-seasinsoi(1);
delta_I=seasanomi-MAanomi;

% Assess linear SST Sensitivity to seasonal insolation
finalfit=polyfit(delta_I,ssti,1);
finalslope=finalfit(1);
SSTSensitivity=finalslope;
SSTsens=[SSTsens,SSTSensitivity];

finalmdl=fitlm(ssti,delta_I,'linear');
SST_SeasInso_R2=finalmdl.Rsquared.Ordinary;
linear_R2=[];
linear_R2=[linear_R2,SST_SeasInso_R2];

[kyear, lat, day]=meshgrid(0:1:150,latmin:0.1:latmax,1:1:365);
MA_Inso=daily_insolation(kyear,lat,day);
MA_Inso=mean(MA_Inso,1);
MA_Inso=mean(MA_Inso,3);
MA_Inso=squeeze(MA_Inso);

% Calculate SST adjustment to convert seasonal to mean annual SST
MA_anom=(MA_Inso(1:151)-MA_Inso(1));
Seas_anom=(Seas_Inso(1:151)-Seas_Inso(1));
SSTadjustment=(Seas_anom'-MA_anom').*(SSTSensitivity);
SSTadjustmenti=interp1(0:150,SSTadjustment,core(:,1),'linear','extrap');

core_MA=core(:,2)-SSTadjustmenti;
end
core_ma2=[core_ma2,core_MA];

end
core_ma=[core(:,1),mean(core_ma2,2),(std(core_ma2,0,2)./sqrt(2))];
SST_sensitivity=[mean(SSTsens),std(SSTsens)];
Start_day=mean(start_day);
End_day=mean(end_day);
Linear_r2=mean(linear_R2);

figure(4)
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

figure(5)
plot(core_ma(:,1),core_ma(:,2),'-o','LineWidth',3)
hold on
plot(core(:,1),core(:,2),'-o','LineWidth',3)
ylabel('SST Anomaly (°C)')
set(gca,'FontSize',12)
xlabel('Age(ka)')
xlim([115 128])
set(gca,'xdir','reverse') 
title('Last Interglacial')
legend('Mean Annual','Original unadjusted data')
end
