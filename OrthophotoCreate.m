close all
clear all
clc
addpath('/Users/cmbaker9/Documents/MTOOLS')
addpath(genpath('/Users/cmbaker9/Documents/Research/swashzone_dynamics/codes/cameras'))
% genpath('/Users/cmbaker9/Documents/Research/Lab_Experiments/codes/CIRN-Quantitative-Coastal-Imaging-Toolbox')
addpath(genpath('/Users/cmbaker9/Documents/Research/Lab_Experiments/codes/CIRN-Quantitative-Coastal-Imaging-Toolbox/X_CoreFunctions'))
projpath = '/Users/cmbaker9/Documents/Research/swashzone_dynamics/';

%% User input

date2pick = '2021-10-13 19'; % format: yyyy-mm-dd HH

%% STEP 1: Load Data and constants

datapath = '/Volumes/i2rgusBob4';
fnames = getfnames(datapath,'');
eTime = str2num(fnames);

mTime = datetime(eTime,'ConvertFrom','epochtime','TicksPerSecond',1000);
mpick = datenum(date2pick,'yyyy-mm-dd HH');
epick = (mpick-datenum(1970,1,1,0,0,0))*(24*3600*1000);

[~,ip] = min(abs(eTime-epick));

% folder
fname=num2str(eTime(ip));

%% Load File Names + Times

idir = [datapath,'/',num2str(eTime(ip))];
cnamesA = getfnames(idir,'CamBob_A');
cnamesB = getfnames(idir,'CamBob_B');
tA = str2num(cnamesA(:,10:end-5));
tdateA=tA/24/3600/1000+datenum(1970,1,1);
tB = str2num(cnamesB(:,10:end-5));
tdateB=tB/24/3600/1000+datenum(1970,1,1);
dt=mode(diff(tA));
ts=tA-tA(1);

%% Load extrinsics and intrinsics

%  Enter the filepath of the saved  IOEO calibration results
ieopath = '/Users/cmbaker9/Documents/Research/swashzone_dynamics/codes/cameras/61SwashI2Rgus/ForChristine2021/';
% EX = readtable([ieopath,'CameraPositions_YawPitchRoll.txt']);
% extrinsics{1} = EX{1,8:13};
% extrinsics{2} = EX{2,8:13};

IN = readtable([ieopath,'Intrinsics-MetashapeFormat/CamBobA.xlsx']);
intrinsics{1} = intrinsicsMetashape2CIRN(IN);

IN = readtable([ieopath,'Intrinsics-MetashapeFormat/CamBobB.xlsx']);
intrinsics{2} = intrinsicsMetashape2CIRN(IN);

% # PhotoID, X, Y, Z, Omega, Phi, Kappa, r11, r12, r13, r21, r22, r23, r31, r32, r33
Rt = importdata([ieopath,'CameraPositions_R.txt']);
[xc(1),yc(1),zc(1),R(:,:,1)] = organizexyzR(Rt.data(1,:));
[xc(2),yc(2),zc(2),R(:,:,2)] = organizexyzR(Rt.data(2,:));
clear Rt IN

%% Convert Local Grid to SPE/SPN

% load('/Volumes/DUNEXinsitu/Example2021Data/1608217202014/STEREOGRIDS/202012171500_grid10cm_1Hz_partial.mat')
dxy = 0.2;
Xtemp = 60:dxy:280;
Ytemp = 535:dxy:710;
X = repmat(Xtemp',[1 length(Ytemp)]);
Y = repmat(Ytemp,[length(Xtemp) 1]);
[ALat, ALon, spN, spE, Y, X] = frfCoord(X, Y);

tide = 0.468;
% Zelev = zeros(3000,4096)+tide;
Zelev = zeros(size(X))+tide;

%% Make Orthos
% Make Directory
% odir = string(fullfile(projpath,'data','cameras'));
% mkdir(string(fullfile(odir,'Datasets',fname,'ORTHOS')))
sdir = '/Volumes/Lacie/DUNEX2021/PROCESSED/CAM/';
sname = datestr(mTime(ip));
sname = regexprep(sname,':','.');
sname = regexprep(sname,'-','.');
sname = regexprep(sname,' ','_');
geomname = ['x',num2str(Xtemp(1)),'-',num2str(Xtemp(end)),'_y',num2str(Ytemp(1)),'-',num2str(Ytemp(end)),'_dxy',num2str(dxy)];
mkdir(string(fullfile(sdir,sname,'ORTHOS',geomname)))
% t = Tiff([camfiles(3).folder,'/CamBob_A_',camfiles(imgrab).name(10:end-5),'.tiff'],'r')
% t = Tiff([camfiles(3).folder,'/CamBob_B_',camfiles(imgrab).name(10:end-7),'80','.tiff'],'r');

%%
oname=['Rectxyz'];
collecttime= datestr(mTime(ip));
save(fullfile(sdir,sname,'ORTHOS',[oname '.mat']),'X','Y','Zelev','collecttime')
display('only reading Bob_A')

% 800 to 2440 to 3799
for k=1:(size(cnamesA,1))
    for ic = 1%1:2
        % camera info
        if ic == 1
            cam = 'Bob_A';
            tstr = num2str(tA(k));
        elseif ic == 2
            cam = 'Bob_B';
            [~,ib] = min(abs(tB-tA(k)));
            tstr = num2str(tB(ib));
        end
        
        % Load Image
        cname=fullfile(idir,['Cam' cam '_' tstr '.tiff']);
        I=imread(cname);
        [Ir]= imageRectificationmod(I,intrinsics{2},[xc(ic) yc(ic) zc(ic)],R(:,:,ic),spE,spN,Zelev,0);
%         figure;
%         subplot(1,2,1); imshow(I);
%         subplot(1,2,2); imshow(Ir);

        % Save Image
        imwrite(Ir,string(fullfile(sdir,sname,'ORTHOS',geomname,['Cam' cam '_' tstr '.tiff'])))
    end
    k
end


%%
 cam = 'Bob_A';
oname = ['Cam' cam '_' sname];
v = VideoWriter(fullfile(sdir,sname,'ORTHOS',geomname,[oname,'R.avi']));
v.FrameRate=2;
open(v)

for k=1:(size(cnamesA,1))
    
   
    tstr = num2str(tA(k));
    
    
    % Load Image
    cname=fullfile(sdir,sname,'ORTHOS',geomname,['Cam' cam '_' tstr '.tiff']);
    I=imread(cname);
    
    
    writeVideo(v,I)

k
end  
   close(v) 
   disp('saved movie')

save(fullfile(sdir,sname,'ORTHOS',geomname,[oname '_xyz.mat']),'X','Y','Zelev','collecttime')
display('only reading Bob_A')

%%
 cam = 'Bob_A';
oname = ['Cam' cam '_' sname];
v = VideoWriter(fullfile(sdir,sname,'ORTHOS',geomname,[oname,'R_short.avi']));
v.FrameRate=2;
open(v)

for k=1:200%(size(cnamesA,1))
    
   
    tstr = num2str(tA(k));
    
    
    % Load Image
    cname=fullfile(sdir,sname,'ORTHOS',geomname,['Cam' cam '_' tstr '.tiff']);
    I=imread(cname);
    
    
    writeVideo(v,I)

k
end  
   close(v) 
   disp('saved movie')

save(fullfile(sdir,sname,'ORTHOS',geomname,[oname '_xyz.mat']),'X','Y','Zelev','collecttime')
display('only reading Bob_A')