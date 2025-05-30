%% Initialize
clear;
close all;
clc
model = 2; % 1 for flange.bin;  2 for MPCB
NoisyCase = 1;

if model == 1
    % Define Geometry
    length = 350;width = 350;height = 50;
    geo=defaultGeometry();                     
    geo.DSD = 816;                              % Distance Source Detector      (mm)
    geo.DSO = 600;                              % Distance Source Origin        (mm)
    % Detector parameters
    geo.nDetector=[512; 512];					% number of pixels              (px)
    geo.dDetector=[0.127; 0.127]; 				% size of each pixel            (mm)
    geo.sDetector=geo.nDetector.*geo.dDetector; % total size of the detector    (mm)
    % Image parameters
    geo.nVoxel=[length;height;width];           % number of voxels              (vx)
    geo.sVoxel=[length*0.127;height*0.127;width*0.127];        % total size of the image       (mm)
    geo.dVoxel=geo.sVoxel./geo.nVoxel;          % size of each voxel            (mm)
    % Load data and generate projections 
    numProjs = 30;
    angles1=linspace(0,2*pi,numProjs);
    angles=[zeros(1,numProjs);angles1;ones(1,numProjs)*pi/180*(45)];
    filename = "flange.bin";
    fid=fopen(filename,'rb');
    eascan = fread(fid, length*width*height, 'float');
    I = single(reshape(eascan, [length,width,height]));
    I = permute(I, [2 3 1]); % 本实验以读入图像转置后作为标准图像，即[x,y,z]->[y,z,x]
    projections=Ax(I,geo,angles,'interpolated');   
    if NoisyCase == 1
        noise_projections=addCTnoise(projections,'Poisson',1e4);
    elseif NoisyCase == 2
        noise_projections=addCTnoise(projections,'Poisson',1e4,'Gaussian',[0,10]);
    else
        noise_projections=(projections);
    end
elseif model == 2 
    % Define Geometry
    length = 256;width = 256;height = 75;
    geo=defaultGeometry();                     
    geo.DSD = 816;                              % Distance Source Detector      (mm)
    geo.DSO = 600;                              % Distance Source Origin        (mm)
    % Detector parameters
    geo.nDetector=[512; 512];					% number of pixels              (px)
    geo.dDetector=[0.127; 0.127]; 				% size of each pixel            (mm)
    geo.sDetector=geo.nDetector.*geo.dDetector; % total size of the detector    (mm)
    % Image parameters
    geo.nVoxel=[length;height;width];           % number of voxels              (vx)
    geo.sVoxel=[length*0.127;height*0.127;width*0.127];        % total size of the image       (mm)
    geo.dVoxel=geo.sVoxel./geo.nVoxel;          % size of each voxel            (mm)
    % Load data and generate projections 
    numProjs = 30;
    angles1=linspace(0,2*pi,numProjs);
    angles=[zeros(1,numProjs);angles1;ones(1,numProjs)*pi/180*(45)];
    filename = "./MPCB.bin";
    fid=fopen(filename,'rb');
    eascan = fread(fid, length*width*height, 'float');
    I = single(reshape(eascan, [length,width,height]));
    I = permute(I, [2 3 1]); % 本实验以读入图像转置后作为标准图像，即[x,y,z]->[y,z,x]
    projections=Ax(I,geo,angles,'interpolated');    
    if NoisyCase == 1
        noise_projections=addCTnoise(projections,'Poisson',1e5);
    else
        noise_projections=(projections);
    end
else
    error(['parameter "' 'model' '" does not exist' ]);
end

%% SART-SR algorithm                               % Symbols in the paper
SART_lambda=0.8;                                   % λ_SART
lambdared=0.9999;                                  % κ_1
maxiter=300;                                       % Max_iter
smooth_normType = [-inf,-inf,-inf,-inf,-0.5];      % [p1,p2,p3,p4,p5]
smooth_lambda = [0.001,0.001,0.001,0.001,0.0012];  % [α1,α2,α3,α4,α5]
u = 0.2;                                           % τ
ng = 4;                                            % TV_iter

qualmeas={'RMSE','CC','MSSIM','UQI'};              % metrics name
[img,errorL2, quality]=POCS_L0_x_y_z(I,noise_projections,geo,angles,maxiter,smooth_lambda,smooth_normType,u,...
                      'TViter',ng,'lambda',SART_lambda,'lambda_red',lambdared,'verbose',1,'QualMeas',qualmeas); % less important.

% save original image, reconstructed image, RSEN, and metrics('RMSE','CC','MSSIM','UQI')
save(['MPCB_SR_',num2str(numProjs),'_',num2str(smooth_normType(1)),'：',num2str(smooth_lambda(1)),'_',num2str(smooth_normType(2)),'：',num2str(smooth_lambda(2)),'_',num2str(smooth_normType(3)),'：',num2str(smooth_lambda(3)),'_',num2str(smooth_normType(4)),'：',num2str(smooth_lambda(4)),...
    '_',num2str(smooth_normType(5)),'：',num2str(smooth_lambda(5)),'_',num2str(ng),'_',num2str(u),'_TAwTV.mat'],'I','img','errorL2','quality');

%% show results
I = permute(I, [1 3 2]);
img = permute(img, [1 3 2]);
N1 = 31;N2 = 128;N3=128; % It's for MPCB
% N1 = 26;N2 = 175;N3=175; % It's for Flange
figure(1),imshow(reshape(img(:, : ,N1),length,width),[0 1]); axis off;%俯视图
figure(2),imshow(reshape(img(:, N2 ,:),length,height)',[0 1]); axis off;%正视图
figure(3),imshow(reshape(img(N3, : ,:),width,height)',[0 1]); axis off;%侧视图

figure
plot(quality(1,:));
title('Evolution of RMSE per iteration')
figure
plot(quality(2,:));title('Evolution of CC per iteration')
figure
plot(quality(3,:));
title('Evolution of MSSIM per iteration')
figure
plot(quality(4,:));
title('Evolution of UQI per iteration')
 
fprintf('SART-SR \n');
disp(['RMSE:', num2str(quality(1,end)),' RSEN:',num2str(errorL2(1,end)) , ', MSSIM:', num2str(quality(3,end)),  ...
     ', CC:', num2str(quality(2,end)),', UQI:', num2str(quality(4,end))]);
