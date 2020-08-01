clc;close all; clear;
%% 

addpath('.\TV'); addpath('.\FCSBL');

% load 'J' (sensitivity matrix), 'V' (measurement matrix), 
% and 'mask' (mask matrix of the round sensing area).
% NOTE: 'J' is normalized in rows;
%       'V' is subtracted with background signal and normalized.
% emt1: high SNR (>30dB); emt3: low SNR (<20dB).
load('.\data\emt3.mat'); 

% load the mask image of the ground truth phantom position.
load('.\data\ref_img.mat');

load('.\data\s11.mat');
%% One-step Inversion with Regularization

% compute difference operators
dr_filt = [1;-1];
dc_filt = [1,-1];
Dr = convmtx2(dr_filt, 64,64) * mask';
Dc = convmtx2(dc_filt, 64,64) * mask';
DtDr = Dr'*Dr;
DtDc = Dc'*Dc;
Lap = (DtDc + DtDr);

% One-step reconstruction
x_onestep_tv = (J'*J +  1e-5* Lap)\J'* V; % TV reg.

% Display Result
DispRecos(mask * ref_img(:), s11, 64,'linear',0);title('Ground Truth Location');

for i = 1:size(x_onestep_tv,2)    
    DispRecos(x_onestep_tv(:,i), s11, 64,'linear',0);title(['TV Reg. f' num2str(i)]);
end

%% Iterative Total Varaition

% set freq (1,2,3,4)
freq = 4; Vs = V(:,freq);

% 'iso'--> isotropic TV
% 'l1' --> l1-based, anisotropic TV
pars.tv = 'iso';
pars.MAXITER = 200;
X_fista_tv = tv_fista(J,Vs,mask,1e-6,-Inf,Inf,pars);

DispRecos(mask * X_fista_tv(:), s11, 64,'linear',0);title(['Iterative TV f' num2str(freq)]);

%% FCSBL

%%%%%%%%%%%%%%%%% Embed the block structure %%%%%%%%%%%%%%%%%
% generate block structure; run only once
%     b_mat = BlkGenerate(Phi, 'type',2,'len',3,'s11',s11);
%     save('.\data\blk.mat','b_mat');

load('.\data\blk.mat');   % load block structure; 

% block embedding, see Eq. (14) in the paper
J = J * b_mat;

% parameter settings
blkStartLoc = 1:9:size(b_mat,2);
learnlambda = 1;
Result =  FCSBL(J, V, blkStartLoc, learnlambda,...
    'LEARNTYPE', 1, 'lambda', 1e-3, 'prune_gamma',1e-3,...
    'MAX_ITERS', 20,'EPSILON', 1e-8, 'PRINT',0);

X_fcsbl = b_mat*Result.x;  
for i = 1:size(X_fcsbl,2)
    DispRecos(X_fcsbl(:,i), s11, 64,'linear',0);title(['FCSBL f' num2str(i)]);
end

% save('.\results\FCSBL3.mat', 'X_fcsbl');