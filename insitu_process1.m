% Code to organize raw in situ data (pressure gages, advs)

close all
clear all
clc
addpath('/Users/cmbaker9/Documents/MTOOLS')

%% INPUT:

% Selected day and hour to select in situ instrumentation
dateutc   = '2021-10-12 18'; % format: yyyy-mm-dd HH

%% STEP 0: Define file structure and conversions

datapath    = '/Volumes/LaCie/DUNEX2021';
presspath   = fullfile([datapath,'_Press'],'Pressure'); % raw pressure files
velpath     = fullfile([datapath,'_Vect'],'Velocities'); % raw vectrino files
PROCpath    = fullfile(datapath,'PROCESSED','STEP1'); % storage location for processed mat files
eval(['!mkdir ',fullfile(datapath,'PROCESSED')])
eval(['!mkdir ',PROCpath])

% conversions
utc2loc     = -datenum(0,0,0,4,0,0); % conversion between utc and local time (for in situ gages)
mbar2psi  	= 0.0145038; % conversion between mbar and psi (for atmospheric pressure)
psi2Pa      = 6894.76; % conversion between psi and Pa
truN2loc    = 71.86;

%% STEP 1: Load survey information

surveypath = '/Users/cmbaker9/Documents/Research/swashzone_dynamics/data/insitu/swash/';
load([surveypath,'insitu_xyz.mat'],'F');
load([surveypath,'insitu_compass.mat'],'Fc');

for k = 1:9
    % pressure gage pip info
    PS(k).pipeE = F(k).lat;
    PS(k).pipeN = F(k).lon;
    PS(k).pipeZ = F(k).elev;
    % adv pip info
    vect(k).pipeE = F(k).lat;
    vect(k).pipeN = F(k).lon;
    vect(k).pipeZ = F(k).elev;
end

%% STEP 1: Pressure gages

% STEP 1a: read and organize pressure data
fnames = dir([presspath,'/DUNEX_2021_Seabird_Pressures*']);
for i = 1:length(fnames)
    ftemp = fnames(i).name;
    tp(i) = datenum(ftemp(:,30:end-4),'yyyy_mm_dd_HHMM');
end

[~,j] = nanmin(abs(tp-datenum(dateutc,'yyyy-mm-dd HH')+utc2loc));
fread = fnames(j).name;

A = importdata(fullfile(presspath,fread),',',4);
tGMT = datenum(A.textdata(5:end,1));

for k = 1:9
    PS(k).V = A.data(:,k+1);
    PS(k).tm = tGMT;
    PS(k).Vrange = 30; % setting
end

% STEP 1b: Read atmospheric data and interpolate to resolution of sensor
barompath   = ['https://chlthredds.erdc.dren.mil/thredds/dodsC/frf/meteorology/'...
        'barometer/derivedBarom/2021/FRF-met_barometer_derivedBarom_',dateutc(1:4),...
        dateutc(6:7),'.nc'];

psensor     = deblank(ncread(barompath,'station_name')');
patmfrf  	= ncread(barompath,'airPressure')*mbar2psi; % raw unit: mbar
ptime       = ncread(barompath,'time')/24/3600+datenum(1970,01,01,0,0,0);
patm        = interp1(ptime,patmfrf,tGMT);

for k = 1:9
    % converted absolute pressure,V= .1+(4.8*P/Vrange)
    PS(k).psia = ((PS(k).V-0.1)*PS(k).Vrange)/4.8;
    % atmospheric conditions
    PS(k).patm = patm; % psi
    PS(k).patmSource = barompath;
    % salinity
    PS(k).salinity = 35; % psu
end

%% STEP 2: Velocities

% STEP 2a: read raw vectrino data
fnames = dir([velpath,'/DUNEX2021_Vectrino_Velocities*']);
for i = 1:length(fnames)
    ftemp = fnames(i).name;
    tp(i) = datenum(ftemp(:,31:end-4),'yyyy_mm_dd_HHMM');
end

[~,j] = nanmin(abs(tp-datenum(dateutc,'yyyy-mm-dd HH')+utc2loc));
fread = fnames(j).name;

A = importdata(fullfile(velpath,fread),',',4);
tGMT = datenum(A.textdata(5:end,1));

for k = 1:9
    % read voltages
    vect(k).xrawV = A.data(:,3*k-1);
    vect(k).yrawV = A.data(:,3*k);
    vect(k).zrawV = A.data(:,3*k+1);
    % convert to velocity
    vect(k).xrawvel = 2.8*(vect(k).xrawV)-7;
    vect(k).yrawvel = 2.8*(vect(k).yrawV)-7;
    vect(k).zrawvel = 2.8*(vect(k).zrawV)-7;
    % store time
    vect(k).tm = tGMT;
    % store compass information
    vect(k).magDec  = 10.9667; % magnetic declination, 10deg 58'
    vect(k).compass = Fc(k).compass;
    vect(k).tilt    = Fc(k).tilt;
    vect(k).roll    = Fc(k).roll;
    loc2NE(k)       = mean(Fc(k).compass)+vect(k).magDec-180;
    loc2FRF(k)      = loc2NE(k)-(90-truN2loc);
    % calculate vect in terms of NE and local FRF XY
    velX            = vect(k).xrawvel;
    velY            = vect(k).yrawvel;
    vect(k).nvel    = velY.*cosd(loc2NE(k))+velX.*sind(loc2NE(k));
    vect(k).evel    = velX.*cosd(loc2NE(k))-velY.*sind(loc2NE(k));
    vect(k).yvel    = velY.*cosd(loc2FRF(k))+velX.*sind(loc2FRF(k));
    vect(k).xvel    = velX.*cosd(loc2FRF(k))-velY.*sind(loc2FRF(k));
    vect(k).heading = loc2NE;
end

%% Save data

eval(['save -v7.3 ',PROCpath,'/PS_',fread(31:end-4),'.mat PS']);
eval(['save -v7.3 ',PROCpath,'/vect_',fread(31:end-4),'.mat vect']);