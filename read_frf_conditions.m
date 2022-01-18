function [hsig,tp,dir,sprd,watlev,f,See,time_waves,time_watlev] = read_frf_conditions(off_buoy,tide_inst)
% This script will read .nc files and provide the measured conditions.
% INPUT:
% off_buoy  = location of .nc file of offshore wave buoy
% tide_inst = location of .nc file of noaa tidal gage
% OUTPUT:
% hsig      = significant wave height measured by offshore wave buoy
% tpfb      = peak period from frequency band? NOTE: currently not outputting this...
% tp        = peak period
% dir       = mean wave direction, aka incidence angle
% sprd      = directional spread
% f         = frequency for spectra
% See       = measured sprectra at wave buoy
% time      = time range for measurments

%% STEP 1: Specify terms of interest

truN2loc  = 71.86; % Offset between true north and the frf coordinates (relative to shore-normal)
param = {'time_waves';'hsig';'tp';'dir';'sprd';'watlev';'f';'See'};
name  = {'time';'waveHs';'waveTp';'waveMeanDirection';...
    'directionalPeakSpread';'waterLevel';'waveFrequency';'waveEnergyDensity'};

%% STEP2: Process 11-m AWAC observations

% STEP2a: load data
for i = 1:length(param(:,1))
    eval([param{i},' = ncread(off_buoy,name{i});'])
end

% STEP2b: manipulation to convert to local time and coordinates
dir    = -((dir)-truN2loc); % convert from true north to FRF coordinates
time_waves   = time_waves/24/3600+datenum(1970,01,01,0,0,0);

% STEP2c: round time to nearest hour
for i = 1:length(time_waves)
    time_waves_round(i,:)  = round(datenum_round_off(time_waves(i),'hour'),9);
end

time_waves = time_waves_round;

%%

% STEP3a: load data
watlev          = ncread(tide_inst,'waterLevel');
time_watlev     = ncread(tide_inst,'time');

% STEP3b: convert values for local time and coordinates
time_watlev     = time_watlev/24/3600+datenum(1970,01,01,0,0,0);

end
