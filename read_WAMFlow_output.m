%% Read WAMFlow output
%  Makes WAMs for i2Rgus Camera Dec2020 Set Up

%% Housekeeping
close all
clear all
addpath('/Users/cmbaker9/Documents/MTOOLS')
addpath(genpath('/Users/cmbaker9/Documents/Research/swashzone_dynamics/codes/cameras'))
genpath('/Users/cmbaker9/Documents/Research/Lab_Experiments/codes/CIRN-Quantitative-Coastal-Imaging-Toolbox')
% addpath(genpath('/Users/cmbaker9/Documents/Research/Lab_Experiments/codes/CIRN-Quantitative-Coastal-Imaging-Toolbox/X_CoreFunctions'))
projpath = '/Users/cmbaker9/Documents/Research/swashzone_dynamics/';

%% User Input

date2pick = '12.Oct.2021_17.59.01'; % format: yyyy-mm-dd HH

% Averaging Window (seconds)
awin=22;

% Overlap (seconds)
ovr=4;

% AVI Flag
aflag=1;

% camera name
cam = 'Bob_A';

%%

odir = '/Volumes/Lacie/DUNEX2021/PROCESSED/CAM/';

dxy = 1;
Xtemp = 60:dxy:280;
Ytemp = 535:dxy:710;

subreg = 1;
if subreg == 1
    Xs = 76:dxy:140;
    Ys = 580:dxy:700;
%     ix = find(Xs==Xtemp);
    [~,ix(1)] = min(abs(Xtemp-Xs(1)));
    [~,ix(2)] = min(abs(Xtemp-Xs(end)));
    [~,iy(1)] = min(abs(Ytemp-Ys(1)));
    [~,iy(2)] = min(abs(Ytemp-Ys(end)));
end

geomname = ['x',num2str(Xtemp(1)),'-',num2str(Xtemp(end)),'_y',num2str(Ytemp(1)),'-',num2str(Ytemp(end)),'_dxy',num2str(dxy)];

sname = date2pick;
idir = fullfile(odir,sname,'ORTHOS',geomname);
cnames = getfnames(idir,['Cam',cam,'_16']);
% change when not mary!
L = cnames(:,10:end-5);
t = str2num(L);

if subreg == 1
    geomname = ['x',num2str(Xs(1)),'-',num2str(Xs(end)),'_y',num2str(Ys(1)),'-',num2str(Ys(end)),'_dxy',num2str(dxy)];
end

oset = ['Cam',cam,'_'];
oname=[num2str(t(1)*1000) '_' oset '_win' num2str(awin) '_ovr' num2str(ovr)];

fname = fullfile(odir,sname,'WAMS',geomname,[oname,'R.22s.4dt.optFlow.']);

%% Read data

for i = 1:416
    A = readtable([fname,'u',num2str(i),'.txt']);
    vel_u(:,:,i) = A{:,:};
    A = readtable([fname,'v',num2str(i),'.txt']);
    vel_v(:,:,i) = A{:,:};
end


%%

dateutc   = '2021-10-12 18'; % format: yyyy-mm-dd HH
[vect,PS] = load_insitu_for_plots(dateutc);

%%

xsel = 120;
ysel = 623;
[~,ix] = nanmin(abs(Xs-xsel));
[~,iy] = nanmin(abs(Ys-ysel));
velu = squeeze(vel_u(iy,ix,:));
velv = squeeze(vel_v(iy,ix,:));
t = 1:3:size(vel_u,3);

%%

wamts = datenum('Oct-12-2021 17:59.01');
wamt = wamts:datenum(0,0,0,0,0,4):wamts+datenum(0,0,0,0,0,4)*415;

figure;
plot(vect(6).tm-datenum(0,0,0,4,0,0),vect(6).xvel,'k')
hold on
plot(vect(6).tm-datenum(0,0,0,4,0,0),vect(6).yvel,'r')
plot(wamt,velu,'k','LineWidth',1.2)
plot(wamt,velv,'r','LineWidth',1.2)

%%
vectxy = movmean(sqrt(vect(6).xvel(:,6).^2+vect(6).yvel(:,6).^2),8);
pl64ta(sqrt(vect(6).xvel(:,6).^2+vect(6).yvel(:,6).^2),16)
velxy = sqrt(velu.^2+velv.^2);
figure;
plot(vect(6).tm-datenum(0,0,0,4,0,0),vectxy,'k')
hold on
plot(wamt,velxy,'r','LineWidth',1.2)