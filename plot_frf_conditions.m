% Script to read FRF data used to identify the 
clear all
close all
clc
addpath('/Users/cmbaker9/Documents/MTOOLS')
genpath = '/Users/cmbaker9/Documents/Research/swashzone_dynamics/';

% FIGURE FOLDER: where to save figures
fsubfolder = datestr(date,'yy-mm-dd');
figfolder = [genpath,'figures/insitu/',fsubfolder,'/'];
eval(['!mkdir ',figfolder])

%% STEP 1: Load Data and constants

subfolder = 'data/insitu/frf/';
datamonth = '202110';
off_buoy  = 'https://chlthredds.erdc.dren.mil/thredds/dodsC/frf/oceanography/waves/awac-11m/2021/FRF-ocean_waves_awac-11m_202110.nc';%[genpath,subfolder,'FRF-ocean_waves_8m-array_',datamonth,'.nc']; % wave buoy at 11 m water depth
tide_inst = 'https://chlthredds.erdc.dren.mil/thredds/dodsC/frf/oceanography/waterlevel/eopNoaaTide/2021/FRF-ocean_waterlevel_eopNoaaTide_202110.nc';%[genpath,subfolder,'FRF-ocean_waterlevel_eopNoaaTide_',datamonth,'.nc']; % noaa gage at pier to use water level

[obs.hsig,obs.tp,obs.dir,obs.sprd,obs.watlev,obs.f,obs.See,obs.time_waves,obs.time_watlev] = read_frf_conditions(off_buoy,tide_inst)

%% Load tidal elevation from NOAA

utc2loc     = -datenum(0,0,0,4,0,0); % conversion between utc and local time (nearshore instruments)

tidepath = [genpath,'data/insitu/frf/CO-OPS_8651370_met.csv'];

opts = delimitedTextImportOptions("NumVariables", 5);

% Specify range and delimiter
opts.DataLines = [2, Inf];
opts.Delimiter = ",";

% Specify column names and types
opts.VariableNames = ["Date", "TimeGMT", "Predictedft", "Preliminaryft", "Verifiedft"];
opts.VariableTypes = ["categorical", "string", "double", "double", "categorical"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Specify variable properties
opts = setvaropts(opts, "TimeGMT", "WhitespaceRule", "preserve");
opts = setvaropts(opts, ["Date", "TimeGMT", "Verifiedft"], "EmptyFieldRule", "auto");

M = readtable(tidepath,opts);

datenum(string(M{:,1}),'yyyy/mm/dd')

noaa.tide = M{:,4};
noaa.time = datenum(append(string(M{:,1}),string(M{:,2})),'yyyy/mm/ddHH:MM');
% for i = 1:size(M,1)
% %     time(i) = datenum(join([datestr(M{i,1}),' ',datestr(M{i,2},'HH:MM:SS.FFF')]));
%     time(i) = datenum(string(M{i,1}),'yyyy/mm/dd') + ,datestr(string(M{i,2}),'HH:MM:SS.FFF')]));
% end

%% List of simulations

% Simulations ran
tpick(1,:) = '13-Oct-2021 20:00:00'; % in UTC
% tpick(2,:) = '07-Oct-2021 14:00:00';
ex.date  = datestr(datenum(tpick)+utc2loc); % sim2: not tried subgrid or exactly right angle

for i = 1:length(ex.date(:,1))
    [~,k] = nanmin(abs(obs.time_waves-datenum(ex.date(i,:))));
    [~,j] = nanmin(abs(obs.time_watlev-datenum(ex.date(i,:))));
    [~,j2] = nanmin(abs(noaa.time-datenum(ex.date(i,:))));
    ex.time(i) = datenum(ex.date(i,:));
    ex.hsig(i) = obs.hsig(k);
    ex.Tp(i)   = obs.tp(k);
    ex.dir(i)  = obs.dir(k);
    ex.sprd(i) = obs.sprd(k);
    ex.watlev(i) = obs.watlev(j);
    ex.noaatide(i) = noaa.tide(j2);
end

%% Figure for presentations
h = [0 0];
xlabsloc = datenum(str2double(datamonth(1:4)),str2double(datamonth(5:6)),1,0,0,0):7:datenum(str2double(datamonth(1:4)),str2double(datamonth(5:6)),29,0,0,0);
xlabs = {[datamonth(5:6),'/01'],[datamonth(5:6),'/08'],[datamonth(5:6),'/15'],[datamonth(5:6),'/22'],[datamonth(5:6),'/29']};

figure('units','inches','position',[1 1 12 10],'Color','w');
subplot(511)
% plot(obs.time_watlev,obs.watlev,'k','LineWidth',1.5)
hold on
plot(noaa.time,noaa.tide,'k','LineWidth',1.5)
hold on
% scatter(ex.time,ex.watlev,100,'filled','MarkerFaceColor','r')
scatter(noaa.time(j2),ex.noaatide,100,'filled','MarkerFaceColor','r')
p1 = plot([obs.time_waves(1) obs.time_waves(end)],h,'-.','Color','k','LineWidth',1.5);
p1.Color(4) = 0.35;
datetick('x','keeplimits')
xlim([obs.time_waves(1) obs.time_waves(end)])
ylim([-1 1.4])
grid on
box on
h1=gca;
set(h1,'tickdir','out','xminortick','on','yminortick','on');
set(h1,'xtick',xlabsloc,'xticklabel',{'','','','',''});
set(h1,'ticklength',1.5*get(h1,'ticklength'));
set(h1,'ydir','normal');
set(h1,'fontsize',16);
ylabel('$\xi~(\mathrm{m})$','interpreter','latex','fontsize',22);
% title('Offshore~Conditions','interpreter','latex','fontsize',22);
% h2=legend('','ex');
% set(h2,'interpreter','latex','fontsize',20,'orientation','horizontal','Location','northeast');

subplot(512)
plot(obs.time_waves,obs.hsig,'k','linewidth',1.5);
hold on
scatter(ex.time,ex.hsig,200,'filled','MarkerFaceColor','r')
xlim([obs.time_waves(1) obs.time_waves(end)])
ylim([0 3.5])
grid on
h1=gca;
set(h1,'tickdir','out','xminortick','on','yminortick','on');
set(h1,'xtick',xlabsloc,'xticklabel',{'','','','','',''});
set(h1,'ticklength',1.5*get(h1,'ticklength'));
set(h1,'ydir','normal'); % may need to % out
set(h1,'fontsize',16);
ylabel('$H_{s}~(\mathrm{m})$','interpreter','latex','fontsize',22);
% h2=legend('8-m','ex');
% set(h2,'interpreter','latex','fontsize',20,'orientation','horizontal','Location','northeast');

subplot(513)
plot(obs.time_waves,obs.tp,'k','linewidth',1.5);
hold on
scatter(ex.time,ex.Tp,200,'filled','MarkerFaceColor','r')
% p1 = plot([obs.time_waves(1) obs.time_waves(end)],h,'-.','Color','k','LineWidth',1.5);
% p1.Color(4) = 0.35;
ylim([0 20]);
xlim([obs.time_waves(1) obs.time_waves(end)])
grid
h1=gca;
set(h1,'tickdir','out','xminortick','on','yminortick','on');
set(h1,'xtick',xlabsloc,'xticklabel',{'','','','','',''})
set(h1,'ticklength',1.5*get(h1,'ticklength'));
set(h1,'ydir','normal'); % may need to % out
set(h1,'fontsize',16);
% xlabel('Date~{(Local)}','interpreter','latex','fontsize',22);
ylabel('$T_p~(s)$','interpreter','latex','fontsize',22);

subplot(514)
plot(obs.time_waves,obs.dir,'k','linewidth',1.5);
hold on
scatter(ex.time,ex.dir,200,'filled','MarkerFaceColor','r')
p1 = plot([obs.time_waves(1) obs.time_waves(end)],h,'-.','Color','k','LineWidth',1.5);
p1.Color(4) = 0.35;
ylim([-40 40]);
xlim([obs.time_waves(1) obs.time_waves(end)])
grid
h1=gca;
set(h1,'tickdir','out','xminortick','on','yminortick','on');
set(h1,'xtick',xlabsloc,'xticklabel',{'','','','','',''})
set(h1,'ticklength',1.5*get(h1,'ticklength'));
set(h1,'ydir','normal'); % may need to % out
set(h1,'fontsize',16);
% xlabel('Date~{(Local)}','interpreter','latex','fontsize',22);
ylabel('$\theta~(^{\circ})$','interpreter','latex','fontsize',22);

subplot(515)
plot(obs.time_waves,obs.sprd,'k','linewidth',1.5);
hold on
scatter(ex.time,ex.sprd,200,'filled','MarkerFaceColor','r')
% p1 = plot([obs.time_waves(1) obs.time_waves(end)],h,'-.','Color','k','LineWidth',1.5);
% p1.Color(4) = 0.35;
ylim([0 40]);
xlim([obs.time_waves(1) obs.time_waves(end)])
grid
h1=gca;
set(h1,'tickdir','out','xminortick','on','yminortick','on');
set(h1,'xtick',xlabsloc,'xticklabel',xlabs)
set(h1,'ticklength',1.5*get(h1,'ticklength'));
set(h1,'ydir','normal'); % may need to % out
set(h1,'fontsize',16);
xlabel('Date~{(Local)}','interpreter','latex','fontsize',22);
ylabel('$\sigma_{\theta}~(^{\circ})$','interpreter','latex','fontsize',22);

sname = [figfolder,'wave_conditions'];%sim',simno,'_eta_vort_PSchanges.pdf'];
print(sname,'-dpng')