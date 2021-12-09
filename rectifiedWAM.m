%% Function Description
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
dxy = 1;
Xtemp = 60:dxy:280;
Ytemp = 535:dxy:710;
geomname = ['x',num2str(Xtemp(1)),'-',num2str(Xtemp(end)),'_y',num2str(Ytemp(1)),'-',num2str(Ytemp(end)),'_dxy',num2str(dxy)];

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
%% Load File Names + Times
odir = '/Volumes/Lacie/DUNEX2021/PROCESSED/CAM/';
sname = date2pick;
idir = fullfile(odir,sname,'ORTHOS',geomname);
cnames = getfnames(idir,['Cam',cam,'_16']);
% change when not mary!
L = cnames(:,10:end-5);
t = str2num(L);
tdate=t/24/3600/1000+datenum(1970,1,1);
dt=mode(diff(t));
ts=(t-t(1))/1000;
oset = ['Cam',cam,'_'];


%% Create Time Stamps and Windows

Tbin= (awin/2):ovr: (ts(end)-awin/2); % Centered On
Tlow=Tbin-awin/2; % Upper Limit
Thigh=Tbin+awin/2; % Lower Limit


% Tbin= (awin/2):ovr: (length(ts)*fr-awin/2); % Centered On
% Tlow=Tbin-awin/2; % Upper Limit
% Thigh=Tbin+awin/2; % Lower Limit

%% Create Structure to Hold Images
Io=imread(fullfile(idir,[oset,L(1,:),'.tiff']));
[r c co]=size(Io);
% Iwam=zeros(r,c,length(Tbin));
if subreg == 0
    Iwam=zeros(c,r,length(Tbin));
elseif subreg == 1
    Iwam=zeros(iy(2)-iy(1)+1,ix(2)-ix(1)+1,length(Tbin));
end
Ncounter=zeros(1,length(Tbin));
clear('Io')


%% Load Images

for k=1:length(L)
    
% Load Image (Gray- Matlab cannot handle RGB array size)
% I=double(rgb2gray(imread(fullfile(idir,[oset,L(k,:),'.tiff']))));
I=double(flipud(rgb2gray(imread(fullfile(idir,[oset,L(k,:),'.tiff'])))'));
if subreg == 1
    I = squeeze(I(iy(1):iy(2),ix(1):ix(2),:));
end

% Find Bins Image Belongs into
tcheck=ts(k);    
lcheck=Tlow-tcheck;
hcheck=Thigh-tcheck;
ind=find(lcheck<=0 & hcheck>0);

% Add Image to Average
Iwam(:,:,ind)=Iwam(:,:,ind)+I;

% Add Counter for average at end
Ncounter(ind)=Ncounter(ind)+1;

k
end

disp('Loaded Images')

%% Take Average
for k=1:length(Tbin)
    Iwam(:,:,k)=Iwam(:,:,k)./Ncounter(k);
end
Iwam=uint8(Iwam);

TbinS=Tbin;
Tbin=Tbin+t(1);
disp('Averaged Images')

if subreg == 1
    geomname = ['x',num2str(Xs(1)),'-',num2str(Xs(end)),'_y',num2str(Ys(1)),'-',num2str(Ys(end)),'_dxy',num2str(dxy)];
end

%% Save File (MAT)
mkdir(string(fullfile(odir,sname,'WAMS',geomname)))

oname=[num2str(t(1)*1000) '_' oset '_win' num2str(awin) '_ovr' num2str(ovr)];
save(fullfile(odir,sname,'WAMS',geomname,[oname 'R.mat']),'Iwam','Ncounter','TbinS','awin','ovr','TbinS')
disp('saved Mat File')

%% Make A Movie
if aflag==1
    
v = VideoWriter(fullfile(odir,sname,'WAMS',geomname,[oname,'R.avi']));
v.FrameRate=2;
open(v)

for k=1:length(Tbin)
writeVideo(v,(Iwam(:,:,k)))

k
end  
   close(v) 
   disp('saved movie')
end
    
    %% make short video
    
    oname=[num2str(t(1)*1000) '_' oset '_win' num2str(awin) '_ovr' num2str(ovr)];
    v = VideoWriter(fullfile(odir,sname,'WAMS',geomname,[oname,'R_short.avi']));
v.FrameRate=2;
open(v)

for k=24:224
writeVideo(v,(Iwam(:,:,k)))

k
end  
   close(v) 
   disp('saved movie')
    
% green and blue
% each color band
    

%% make short video
geomname = ['x',num2str(Xtemp(1)),'-',num2str(Xtemp(end)),'_y',num2str(Ytemp(1)),'-',num2str(Ytemp(end)),'_dxy',num2str(dxy)];

oname=[num2str(t(1)*1000) '_' oset '_win' num2str(awin) '_ovr' num2str(ovr)];
    v = VideoWriter(fullfile(odir,sname,'ORTHOS',geomname,[oname,'R_short.avi']));
    v.FrameRate=8;
    open(v)
    
    for k=24:124
        I=imread(fullfile(idir,[oset,L(k,:),'.tiff']));
        writeVideo(v,I)
        
        k
    end
    close(v)
    disp('saved movie')
    
    % green and blue
    % each color band
    

