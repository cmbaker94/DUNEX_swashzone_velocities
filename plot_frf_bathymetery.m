% Read and plot the FRF bathymetry
clc
clear all
close all
addpath(genpath('/Users/cmbaker9/Documents/MTOOLS')) % location where cmocean is saved

%% STEP 0: User input:

bathydate   = '20210928'; % bathymetry date
clim        = [-8.5 8.5]; % color map limits (m)
cmres       = 40; % color map resolution
conres      = 1; % contour resolution (m)

%% STEP 1: Load Data

% website with bathymetry
bathy   = ['https://chlthredds.erdc.dren.mil/thredds/dodsC/frf/geomorphology/',...
    'DEMs/surveyDEM/data/FRF_geomorphology_DEMs_surveyDEM_',bathydate,'.nc'];

x       = ncread(bathy,'xFRF'); % cross-shore
y       = ncread(bathy,'yFRF'); % alongshore
h       = ncread(bathy,'elevation')'; % elevaiton

%% transects

transects = ['https://chlthredds.erdc.dren.mil/thredds/dodsC/frf/geomorphology/elevationTransects/survey/data/FRF_geomorphology_elevationTransects_survey_20211021.nc']
x       = ncread(transects,'xFRF'); % cross-shore
y       = ncread(transects,'yFRF'); % cross-shore
elev    = ncread(transects,'elevation'); % cross-shore

%% STEP 2: Create plot

figure('units','inches','position',[1 1 8 11],'Color','w');
axes('Position',[0.15 0.12 0.75 0.75])
contourf(x,y,h,200,'linestyle','none');
hold on
for i = 1:conres:abs(clim(1))
    contour(x,y,-h,[i,i],'Color',[.5 .5 .5],'LineWidth',1)
end
text(430,1265,'$z$ (m)','interpreter','latex','fontsize',23);
xlim([min(x,[],'all') max(x,[],'all')])
ylim([min(y,[],'all') max(y,[],'all')])
ylabel('$y$ (m)','interpreter','latex','fontsize',23);
xlabel('$x$ (m)','interpreter','latex','fontsize',23);
colormap(cmocean('topo',cmres))
caxis(clim)
cb = colorbar('Position', [0.15 0.89 0.75 0.025],'Location','north');
h1=gca;
set(h1,'tickdir','out','xminortick','on','yminortick','on');
set(h1,'ticklength',1.5*get(h1,'ticklength'));
set(h1,'ydir','normal');
set(h1,'fontsize',20);

Sname = ['FRF_bathymetry_',bathydate];
print(Sname,'-dpng')
