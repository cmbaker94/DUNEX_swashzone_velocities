% Code to convert pressures gages and plot measurments

close all
clear all
clc
addpath('/Users/cmbaker9/Documents/MTOOLS')
genpath = '/Users/cmbaker9/Documents/Research/swashzone_dynamics/';

% FIGURE FOLDER: where to save figures
fsubfolder = datestr(date,'yy-mm-dd');
figfolder = [genpath,'figures/insitu/',fsubfolder,'/'];
eval(['!mkdir ',figfolder])

%% INPUT:

% Selected day and hour to select in situ instrumentation
dateutc   = '2021-10-12 18'; % format: yyyy-mm-dd HH
utc2loc   = -datenum(0,0,0,4,0,0); % conversion between utc and local time (for in situ gages)
psi2Pa      = 6894.76; % conversion between psi and Pa

%% STEP 0: Define file structure and conversions

datapath    = '/Volumes/LaCie/DUNEX2021';
PROCpath1    = fullfile(datapath,'PROCESSED','STEP1'); % storage location for processed mat files
PROCpath2    = fullfile(datapath,'PROCESSED','STEP2'); % storage location for processed mat files

%% STEP 1: Load data from step1 of processing

% STEP 1a: read and organize pressure data
fnames = dir([PROCpath1,'/PS*']);
for i = 1:length(fnames)
    ftemp = fnames(i).name;
    tp(i) = datenum(ftemp(4:end-4),'yyyy_mm_dd_HHMM');
end

[~,j] = nanmin(abs(tp-datenum(dateutc,'yyyy-mm-dd HH')+utc2loc));
fread = fnames(j).name;

load([PROCpath1,'/PS_',fread(4:end-4),'.mat'])
load([PROCpath1,'/vect_',fread(4:end-4),'.mat'])

%% STEP 2: reconstruct water surface

% Convert psia to depth, assuming hydrostatic
g = 9.81; % gravity, m/s
for k = 1:9
    if PS(k).salinity == 35
        rho = 1026; % kg/m3
    end
    press = (PS(k).psia-PS(k).patm)*psi2Pa; % convert to metric, Pa (or N/m2)
    pressmedfilt = medfilt1(press,5);
    h(k,:) = pressmedfilt/(rho*g); 
    eta(k,:) = detrend(h(k,:));
end



%% Plot

xprange = [PS(2).tm(1)+datenum(0,0,0,0,5,0) PS(2).tm(1)+datenum(0,0,0,0,10,0)];

figure('units','inches','position',[1 1 20 12],'Color','w');
subplot(331)
plot(PS(1).tm,h(1,:),'k','LineWidth',0.3)
hold on
xlim(xprange)
datetick('x','keeplimits')
ylim([-0.1 -0.03])
grid on
h1=gca;
set(h1,'tickdir','out','xminortick','on','yminortick','on');
% set(h1,'xtick',xlabsloc,'xticklabel',{'','','','',''});
set(h1,'ticklength',1.5*get(h1,'ticklength'));
set(h1,'fontsize',16);
ylabel('$h~(\mathrm{m})$','interpreter','latex','fontsize',16);
title('\bf{F1}','interpreter','latex','fontsize',16);

subplot(332)
plot(PS(2).tm,h(2,:),'k','LineWidth',0.3)
hold on
xlim(xprange)
datetick('x','keeplimits')
ylim([-0.1 0.1])
grid on
h1=gca;
set(h1,'tickdir','out','xminortick','on','yminortick','on');
% set(h1,'xtick',xlabsloc,'xticklabel',{'','','','',''});
set(h1,'ticklength',1.5*get(h1,'ticklength'));
set(h1,'fontsize',16);
title('\bf{F2}','interpreter','latex','fontsize',16);

subplot(333)
plot(PS(3).tm,h(3,:),'k','LineWidth',0.3)
hold on
xlim(xprange)
datetick('x','keeplimits')
ylim([-0.2 1])
grid on
h1=gca;
set(h1,'tickdir','out','xminortick','on','yminortick','on');
% set(h1,'xtick',xlabsloc,'xticklabel',{'','','','',''});
set(h1,'ticklength',1.5*get(h1,'ticklength'));
set(h1,'fontsize',16);
title('\bf{F3}','interpreter','latex','fontsize',16);

subplot(334)
plot(PS(4).tm,h(4,:),'k','LineWidth',0.3)
hold on
xlim(xprange)
datetick('x','keeplimits')
ylim([-0.1 -0.03])
grid on
h1=gca;
set(h1,'tickdir','out','xminortick','on','yminortick','on');
% set(h1,'xtick',xlabsloc,'xticklabel',{'','','','',''});
set(h1,'ticklength',1.5*get(h1,'ticklength'));
set(h1,'fontsize',16);
ylabel('$h~(\mathrm{m})$','interpreter','latex','fontsize',16);
title('\bf{F4}','interpreter','latex','fontsize',16);

subplot(335)
plot(PS(5).tm,h(5,:),'k','LineWidth',0.3)
hold on
xlim(xprange)
datetick('x','keeplimits')
ylim([-0.1 0.1])
grid on
h1=gca;
set(h1,'tickdir','out','xminortick','on','yminortick','on');
% set(h1,'xtick',xlabsloc,'xticklabel',{'','','','',''});
set(h1,'ticklength',1.5*get(h1,'ticklength'));
set(h1,'fontsize',16);
title('\bf{F5}','interpreter','latex','fontsize',16);

subplot(336)
plot(PS(6).tm,h(6,:),'k','LineWidth',0.3)
hold on
xlim([PS(6).tm(1) PS(6).tm(end)])
datetick('x','keeplimits')
ylim([-0.2 1])
grid on
h1=gca;
set(h1,'tickdir','out','xminortick','on','yminortick','on');
% set(h1,'xtick',xlabsloc,'xticklabel',{'','','','',''});
set(h1,'ticklength',1.5*get(h1,'ticklength'));
set(h1,'fontsize',16);
title('\bf{F6}','interpreter','latex','fontsize',16);

subplot(337)
plot(PS(7).tm,h(7,:),'k','LineWidth',0.3)
hold on
xlim(xprange)
datetick('x','keeplimits')
ylim([-0.1 -0.03])
grid on
h1=gca;
set(h1,'tickdir','out','xminortick','on','yminortick','on');
% set(h1,'xtick',xlabsloc,'xticklabel',{'','','','',''});
set(h1,'ticklength',1.5*get(h1,'ticklength'));
set(h1,'fontsize',16);
ylabel('$h~(\mathrm{m})$','interpreter','latex','fontsize',16);
xlabel('$t~(\mathrm{MM:HH})$','interpreter','latex','fontsize',16);
title('\bf{F7}','interpreter','latex','fontsize',16);

subplot(338)
plot(PS(8).tm,h(8,:),'k','LineWidth',0.3)
hold on
xlim(xprange)
datetick('x','keeplimits')
ylim([-0.1 0.1])
grid on
h1=gca;
set(h1,'tickdir','out','xminortick','on','yminortick','on');
% set(h1,'xtick',xlabsloc,'xticklabel',{'','','','',''});
set(h1,'ticklength',1.5*get(h1,'ticklength'));
set(h1,'fontsize',16);
xlabel('$t~(\mathrm{MM:HH})$','interpreter','latex','fontsize',16);
title('\bf{F8}','interpreter','latex','fontsize',16);

subplot(339)
plot(PS(9).tm,h(9,:),'k','LineWidth',0.3)
hold on
xlim(xprange)
datetick('x','keeplimits')
ylim([-0.2 1])
grid on
h1=gca;
set(h1,'tickdir','out','xminortick','on','yminortick','on');
% set(h1,'xtick',xlabsloc,'xticklabel',{'','','','',''});
set(h1,'ticklength',1.5*get(h1,'ticklength'));
set(h1,'fontsize',16);
xlabel('$t~(\mathrm{MM:HH})$','interpreter','latex','fontsize',16);
title('\bf{F9}','interpreter','latex','fontsize',16);

sname = [figfolder,'PS_',dateutc(1:10),'_',dateutc(12:13)];%sim',simno,'_eta_vort_PSchanges.pdf'];
print(sname,'-dpng')

%%

figure('units','inches','position',[1 1 20 12],'Color','w');
subplot(331)
plot(vect(1).tm,vect(1).xvel,'k','LineWidth',0.3)
hold on
plot(vect(1).tm,vect(1).yvel,'r','LineWidth',0.3)
xlim(xprange)
ylim([-8 8])
datetick('x','keeplimits')
grid on
h1=gca;
set(h1,'tickdir','out','xminortick','on','yminortick','on');
% set(h1,'xtick',xlabsloc,'xticklabel',{'','','','',''});
set(h1,'ticklength',1.5*get(h1,'ticklength'));
set(h1,'fontsize',16);
ylabel('$u~(\mathrm{m/s})$','interpreter','latex','fontsize',16);
title('\bf{F1}','interpreter','latex','fontsize',16);
h2=legend('$u_x$','$u_y$');
set(h2,'interpreter','latex','fontsize',20,'orientation','horizontal','Location','northeast');


subplot(332)
plot(vect(2).tm,vect(2).xvel,'k','LineWidth',0.3)
hold on
plot(vect(2).tm,vect(2).yvel,'r','LineWidth',0.3)
xlim(xprange)
datetick('x','keeplimits')
ylim([-8 8])
grid on
h1=gca;
set(h1,'tickdir','out','xminortick','on','yminortick','on');
% set(h1,'xtick',xlabsloc,'xticklabel',{'','','','',''});
set(h1,'ticklength',1.5*get(h1,'ticklength'));
set(h1,'fontsize',16);
title('\bf{F2}','interpreter','latex','fontsize',16);

subplot(333)
plot(vect(3).tm,vect(3).xvel,'k','LineWidth',0.3)
hold on
plot(vect(3).tm,vect(3).yvel,'r','LineWidth',0.3)
xlim(xprange)
datetick('x','keeplimits')
ylim([-8 8])
grid on
h1=gca;
set(h1,'tickdir','out','xminortick','on','yminortick','on');
% set(h1,'xtick',xlabsloc,'xticklabel',{'','','','',''});
set(h1,'ticklength',1.5*get(h1,'ticklength'));
set(h1,'fontsize',16);
title('\bf{F3}','interpreter','latex','fontsize',16);

subplot(334)
plot(vect(4).tm,vect(4).xvel,'k','LineWidth',0.3)
hold on
plot(vect(4).tm,vect(4).yvel,'r','LineWidth',0.3)
xlim(xprange)
datetick('x','keeplimits')
% ylim([-0.082 -0.03])
grid on
h1=gca;
set(h1,'tickdir','out','xminortick','on','yminortick','on');
% set(h1,'xtick',xlabsloc,'xticklabel',{'','','','',''});
set(h1,'ticklength',1.5*get(h1,'ticklength'));
set(h1,'fontsize',16);
ylabel('$u~(\mathrm{m/s})$','interpreter','latex','fontsize',16);
title('\bf{F4}','interpreter','latex','fontsize',16);

subplot(335)
plot(vect(5).tm,vect(5).xvel,'k','LineWidth',0.3)
hold on
plot(vect(5).tm,vect(5).yvel,'r','LineWidth',0.3)
xlim(xprange)
datetick('x','keeplimits')
ylim([-8 8])
grid on
h1=gca;
set(h1,'tickdir','out','xminortick','on','yminortick','on');
% set(h1,'xtick',xlabsloc,'xticklabel',{'','','','',''});
set(h1,'ticklength',1.5*get(h1,'ticklength'));
set(h1,'fontsize',16);
title('\bf{F5}','interpreter','latex','fontsize',16);

subplot(336)
plot(vect(6).tm,vect(6).xvel,'k','LineWidth',0.3)
hold on
plot(vect(6).tm,vect(6).yvel,'r','LineWidth',0.3)
xlim(xprange)
datetick('x','keeplimits')
ylim([-8 8])
grid on
h1=gca;
set(h1,'tickdir','out','xminortick','on','yminortick','on');
% set(h1,'xtick',xlabsloc,'xticklabel',{'','','','',''});
set(h1,'ticklength',1.5*get(h1,'ticklength'));
set(h1,'fontsize',16);
title('\bf{F6}','interpreter','latex','fontsize',16);

subplot(337)
plot(vect(7).tm,vect(7).xvel,'k','LineWidth',0.3)
hold on
plot(vect(7).tm,vect(7).yvel,'r','LineWidth',0.3)
xlim(xprange)
datetick('x','keeplimits')
ylim([-8 8])
grid on
h1=gca;
set(h1,'tickdir','out','xminortick','on','yminortick','on');
% set(h1,'xtick',xlabsloc,'xticklabel',{'','','','',''});
set(h1,'ticklength',1.5*get(h1,'ticklength'));
set(h1,'fontsize',16);
ylabel('$u~(\mathrm{m/s})$','interpreter','latex','fontsize',16);
xlabel('$t~(\mathrm{MM:HH})$','interpreter','latex','fontsize',16);
title('\bf{F7}','interpreter','latex','fontsize',16);

subplot(338)
plot(vect(8).tm,vect(8).xvel,'k','LineWidth',0.3)
hold on
plot(vect(8).tm,vect(8).yvel,'r','LineWidth',0.3)
xlim(xprange)
datetick('x','keeplimits')
ylim([-8 8])
grid on
h1=gca;
set(h1,'tickdir','out','xminortick','on','yminortick','on');
% set(h1,'xtick',xlabsloc,'xticklabel',{'','','','',''});
set(h1,'ticklength',1.5*get(h1,'ticklength'));
set(h1,'fontsize',16);
xlabel('$t~(\mathrm{MM:HH})$','interpreter','latex','fontsize',16);
title('\bf{F8}','interpreter','latex','fontsize',16);

subplot(339)
plot(vect(9).tm,vect(9).xvel,'k','LineWidth',0.3)
hold on
plot(vect(9).tm,vect(9).yvel,'r','LineWidth',0.3)
xlim(xprange)
datetick('x','keeplimits')
ylim([-8 8])
grid on
h1=gca;
set(h1,'tickdir','out','xminortick','on','yminortick','on');
% set(h1,'xtick',xlabsloc,'xticklabel',{'','','','',''});
set(h1,'ticklength',1.5*get(h1,'ticklength'));
set(h1,'fontsize',16);
xlabel('$t~(\mathrm{MM:HH})$','interpreter','latex','fontsize',16);
title('\bf{F9}','interpreter','latex','fontsize',16);

sname = [figfolder,'vect_',dateutc(1:10),'_',dateutc(12:13)];%sim',simno,'_eta_vort_PSchanges.pdf'];
print(sname,'-dpng')