% This code will read the csv file of gps data, save the relevant
% measurments, and create a few plots of the compass measurments.

clear all
close all
clc
addpath(genpath('/Users/cmbaker9/Documents/MTOOLS'))
addpath('/Users/cmbaker9/Documents/Research/swashzone_dynamics/codes/cameras/61SwashI2Rgus/Toolbox')

%% Load compass observations

genpath = '/Volumes/LaCie/DUNEX2021/RAW/INSITU/GPS/';
fname = [genpath,'GPSCPMAPSS_20210906_1730_ALL_1sec.csv'];
[c1s] = read_compass_data(fname);

fname = [genpath,'GPSCPMAPSS_20210906_1730_ALL_ALL.csv'];
[craw] = read_compass_data(fname);

% convert to local time
craw.time = craw.time+datenum(0,0,0,4,0,0);

%% Load location of in situ gages

surveypath = '/Users/cmbaker9/Documents/Research/swashzone_dynamics/data/insitu/swash/';
load([surveypath,'insitu_xyz.mat'],'F');

for i = 1:9
    Flat(i) = F(i).lat;
    Flon(i) = F(i).lon;
end

%% Define in situ locations: interpretation when compass is 'steady'

iCmin = c1s.compass>193;
iCmax = c1s.compass<198;
iC = (iCmin+iCmax)==2;

inst(1).icompass = [12:29];
inst(2).icompass = [53:66];
inst(3).icompass = [87:99];
inst(4).icompass = [131:140];
inst(5).icompass = [166:175];
inst(6).icompass = [196:204];
inst(7).icompass = [235:242];
inst(8).icompass = [259:266];
inst(9).icompass = [291:297];

for i = 1:9
    Fc(i).compass = c1s.compass(inst(i).icompass);
    Fc(i).tilt = c1s.tilt(inst(i).icompass);
    Fc(i).roll = c1s.roll(inst(i).icompass);
end

saveloc = '/Users/cmbaker9/Documents/Research/swashzone_dynamics/data/insitu/swash/';
eval(['save -v7.3 ',saveloc,'insitu_compass.mat Fc']);

%% Plot lat/lon from compass and known location of in situ gages

figure('units','inches','position',[1 1 8 8],'Color','w');
scatter(c1s.lon,c1s.lat,60,c1s.time,'fill')
hold on
scatter(craw.lon,craw.lat,10,craw.time,'fill')
hold on
scatter(Flon,Flat,100,'k','MarkerFaceColor','r')
xlabel('Lon ($^{\circ}$)','interpreter','latex','fontsize',23);
ylabel('Lat ($^{\circ}$)','interpreter','latex','fontsize',23);
axis equal
box on
grid on
colormap(cmocean('ice'))
labels = {datestr(c1s.time(1)), datestr(c1s.time(length(c1s.time)/2)), datestr(c1s.time(end))};
cb = colorbar('Ticks',[c1s.time(1) c1s.time(length(c1s.time)/2) c1s.time(end)],'TickLabels',labels)
h1=gca;
set(h1,'tickdir','out','xminortick','on','yminortick','on');
set(h1,'ticklength',1.5*get(h1,'ticklength'));
set(h1,'ydir','normal');
set(h1,'fontsize',16);
h2 = legend('Compass Readings','In Situ Gages');
set(h2,'interpreter','latex','fontsize',18);
set(h2,'orientation','vertical','Location','northwest')

sfolder = '/Users/cmbaker9/Documents/Research/swashzone_dynamics/figures/insitu/';
Sname = [sfolder,'compass_readings'];
print(Sname,'-dpng')

%% Plot timeseries of compass measurements

h = [0 0];
xlabs = [datestr(c1s.time(1)); datestr(c1s.time(round(length(c1s.time)/3))); datestr(c1s.time(round(2*length(c1s.time)/3))); datestr(c1s.time(end))];
xlabs = xlabs(:,13:end);
xlabsloc = [c1s.time(1) c1s.time(round(length(c1s.time)/3)) c1s.time(round(2*length(c1s.time)/3)) c1s.time(end)];

figure('units','inches','position',[1 1 12 10],'Color','w');
subplot(311)
scatter(craw.time,craw.compass,10,'MarkerEdgeColor',[0.5 0.6 0.8],'LineWidth',0.5)
hold on
scatter(c1s.time,c1s.compass,20,'MarkerEdgeColor','k','LineWidth',1.5)
p1 = plot([c1s.time(1) c1s.time(end)],h,'-.','Color','k','LineWidth',1.5);
p1.Color(4) = 0.35;
datetick('x','keeplimits')
xlim([c1s.time(1) c1s.time(end)])
box on
ylim([100 360])
grid
h1=gca;
set(h1,'tickdir','out','xminortick','on','yminortick','on');
set(h1,'xtick',xlabsloc,'xticklabel',{'','','','',''});
set(h1,'ticklength',1.5*get(h1,'ticklength'));
set(h1,'ydir','normal');
set(h1,'fontsize',16);
ylabel('Compass ($^{\circ}$)','interpreter','latex','fontsize',22);

subplot(312)
scatter(craw.time,craw.tilt,10,'MarkerEdgeColor',[0.5 0.6 0.8],'LineWidth',0.5)
hold on
scatter(c1s.time,c1s.tilt,20,'MarkerEdgeColor','k','LineWidth',1.5)
p1 = plot([c1s.time(1) c1s.time(end)],h,'-.','Color','k','LineWidth',1.5);
p1.Color(4) = 0.35;
datetick('x','keeplimits')
xlim([c1s.time(1) c1s.time(end)])
ylim([-100 100])
grid on 
box on
h1=gca;
set(h1,'tickdir','out','xminortick','on','yminortick','on');
set(h1,'xtick',xlabsloc,'xticklabel',{'','','','',''});
set(h1,'ticklength',1.5*get(h1,'ticklength'));
set(h1,'ydir','normal');
set(h1,'fontsize',16);
ylabel('Tilt ($^{\circ}$)','interpreter','latex','fontsize',22);

subplot(313)
scatter(craw.time,craw.roll,10,'MarkerEdgeColor',[0.5 0.6 0.8],'LineWidth',0.5)
hold on
scatter(c1s.time,c1s.roll,20,'MarkerEdgeColor','k','LineWidth',1.5)
p1 = plot([c1s.time(1) c1s.time(end)],h,'-.','Color','k','LineWidth',1.5);
p1.Color(4) = 0.35;
ylim([-200 200]);
xlim([c1s.time(1) c1s.time(end)])
grid on
box on
h1=gca;
set(h1,'tickdir','out','xminortick','on','yminortick','on');
set(h1,'xtick',xlabsloc,'xticklabel',xlabs)
set(h1,'ticklength',1.5*get(h1,'ticklength'));
set(h1,'ydir','normal'); % may need to % out
set(h1,'fontsize',16);
xlabel('Date~{(Local)}','interpreter','latex','fontsize',22);
ylabel('Roll ($^{\circ}$)','interpreter','latex','fontsize',22);

sfolder = '/Users/cmbaker9/Documents/Research/swashzone_dynamics/figures/insitu/';
Sname = [sfolder,'compass_tilt_roll_deg'];
print(Sname,'-dpng')

%% Plot timeseries of just the 'compass' data

figure('units','inches','position',[1 1 18 8],'Color','w');
plot(craw.time,craw.compass)
hold on
p1 = plot([c1s.time(1) c1s.time(end)],h,'-.','Color','k','LineWidth',1.5);
p1.Color(4) = 0.35;
ylim([110 350]);
xlim([c1s.time(1) c1s.time(end)])
grid on
box on
h1=gca;
set(h1,'tickdir','out','xminortick','on','yminortick','on');
set(h1,'xtick',xlabsloc,'xticklabel',xlabs)
set(h1,'ticklength',1.5*get(h1,'ticklength'));
set(h1,'ydir','normal'); % may need to % out
set(h1,'fontsize',16);
xlabel('Date~{(Local)}','interpreter','latex','fontsize',22);
ylabel('Compass ($^{\circ}$)','interpreter','latex','fontsize',22);

sfolder = '/Users/cmbaker9/Documents/Research/swashzone_dynamics/figures/insitu/';
Sname = [sfolder,'compass_deg'];
print(Sname,'-dpng')