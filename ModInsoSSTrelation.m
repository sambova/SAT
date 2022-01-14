%% Modern insolation-SST relationships
%% About
%
% This MATLAB code assesses the modern insolation-sea surface temeprature
% relationship in the WPWP, but could be modified to assess the same
% relationship in any region. As written the code reproduces panel b in
% Figure 1 of Bova et al. (2022).
%
% SST data referenced here is the NOAA daily optimum interpolation sea
% surface temperature dataset available at : 
% https://www.ncei.noaa.gov/products/optimum-interpolation-sst
% 
% The matlab script daily_insolation.m to compute daily average insolation 
% written by Huybers and Eisenman (2006) is required to run this code 
% (http://eisenman.ucsd.edu/code.html).
%
% Written by Samantha Bova (sbova@sdsu.edu), Dec, 2021. 
%
% CITATIONS: Bova, S., Rosenthal, Y., Liu, Z., Yan, M., Broccoli, A.J., 
% Godad, S.P., Zeng, C., and Zheng W. Reply to Matters Arising: SAT method
% precludes the reconstruction of interglacial thermal maxima. Nature 
% Matters Arising (2022).
%
% Bova, S., Rosenthal, Y., Liu, Z. et al. Seasonal origin of 
% the thermal maxima at the Holocene and the last interglacial. Nature 589, 
% 548?553 (2021). https://doi.org/10.1038/s41586-020-03155-x
%% Figure 1, panel b 
ncdisp('sst.day.mean.ltm.1971-2000.nc')% W. Pacific 'x73.165.242.58.110.8.44.2.nc'
%%
time = ncread('sst.day.mean.ltm.1971-2000.nc','time'); % s. atlantic x73.165.242.58.124.18.10.10.nc
lat= ncread('sst.day.mean.ltm.1971-2000.nc','lat');%S. pacificx140.172.38.222.126.7.20.19.nc
lon = ncread('sst.day.mean.ltm.1971-2000.nc','lon');%SH  0_180 x73.165.242.58.126.8.57.25.nc
sst = ncread('sst.day.mean.ltm.1971-2000.nc','sst');%Southern'x73.165.242.58.126.8.57.25.nc'
%% Extract SST for WPWP region
%sst = long, lat, time
data=sst(569:640,344:377,:); %-4.125 to 4.125 N, 142.1250 to 159.8750
%% remove no data (land)
[r,c,v] = ind2sub(size(data),find(data == -9.969209968386869e+36));
data(r,c,v)=NaN;
%% spatial average
data=mean(data,1,'omitnan');
data=mean(data,2,'omitnan');
data=squeeze(data);
%% daily insolation using script from Ian Eisenman and Peter Huybers
latitude=0;
latmin=latitude-1;
latmax=latitude+1;
[kyear lat2 day]=meshgrid(0,latmin:0.1:latmax,1:1:365);
Fsw=daily_insolation(kyear,lat2,day);
Fsw=squeeze(Fsw);
Fsw=mean(Fsw,1);
%% Determine the lag between insolation and SST response
% lag tested up to 90 days 
inso2yr=[Fsw,Fsw];
sst_2year=[data;data];
day=1:730;
insofit=[];
for a=1:90;
    inso2yrshift=[inso2yr(end-a:end),inso2yr(1:end-(a+1))];
    temp=corr(inso2yrshift(1:365)',movmean(sst_2year(1:365),20));
    insofit=[insofit,temp];
end
ind=max(insofit)
ind=find(insofit==ind)
a=ind;
inso2yrshift=[inso2yr(end-a:end),inso2yr(1:end-(a+1))];
%% Plot insolation shifted by determined lag with SST
close all
figure(4)
plot(day+a,inso2yr,'-','LineWidth',3)
ylabel('Insolation (W/m2)')
yyaxis right
plot(day,movmean(sst_2year,20),'-','LineWidth',3)
ylabel('SST (°C)')
xlabel ('Days in 2 years')

figure(5)
plot([inso2yrshift(1:365)],movmean(sst_2year(1:365),20),'o')
xlabel('Insolation')
ylabel('SST')
ylim([28 33])
%% fit linear regression to 1 yr inso (shifted by lag) and sst 
WPWPinso=inso2yrshift(1:365)';
WPWPsst=movmean(sst_2year(1:365),20);
p=polyfit(WPWPinso,WPWPsst,1)
yfit=polyval(p,WPWPinso);

close all
figure(1)
plot(WPWPinso,WPWPsst,'.k')
hold on 
plot(WPWPinso,yfit,'LineWidth',3)
ylim([28 33])
ylabel('SST (°C)')
xlabel('Insolation (W/m2)')

yresid = WPWPsst - yfit;
SSresid = sum(yresid.^2);
SStotal = (length(WPWPsst)-1) * var(WPWPsst);
rsq = 1 - SSresid/SStotal
rsq_adj = 1 - SSresid/SStotal * (length(WPWPsst)-1)/(length(WPWPsst)-length(p))
maxdevlin=max(yresid)
