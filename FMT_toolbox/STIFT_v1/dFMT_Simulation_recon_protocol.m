%% Instruction:
% a new STIFT protocol

%% clear memory
clearvars;
close all; 

%% Configurator
exp_dFMT = STIFT;
exp_dFMT.machine_mode = 1;% machine mode
exp_dFMT.configurationSTIFT;

%% Simulator
exp_dFMT.calcuSystemKSTIFT;
exp_dFMT.exRunSTIFT(exp_dFMT.system_K);
exp_dFMT.emRunSTIFT(exp_dFMT.system_K); 

%% System Builder
exp_dFMT.geneMeaMaskSTIFT;
exp_dFMT.geneSolMaskSTIFT;
%% Weighting Matrix
exp_dFMT.geneWMSTIFT(exp_dFMT.system_K);

save([exp_dFMT.data_buffer_directory '/exp_dFMT.mat']);

%% Solver

load('test/exp_dFMT.mat');
%%

close all;
exp_dFMT.inversSTIFT;
exp_dFMT.drawExpSettingSTIFT;
    

exp_dFMT.drawExDetSTIFT;
exp_dFMT.drawEmDetSTIFT;
exp_dFMT.drawFluoSettingSTIFT;
exp_dFMT.drawReconSTIFT;




%% show ground truth (code has error)
%             
% gd =  load("ground_truth.mat");
% iso_int = 0.4;
% ground_truth = gd.fluo_g_disp;
% obj_bound = max(exp_dFMT.node);
% x_display_T = permute(ground_truth,[2,1,3]);
% [m,n,p] = size(x_display_T);
% dl = [obj_bound(2)  obj_bound(1) obj_bound(3)]./[m-1 n-1 p-1];
% [X,Y,Z] = meshgrid(0:(n-1),0:(m-1),0:(p-1));
% X=X*dl(2);
% Y=Y*dl(1);
% Z=Z*dl(3);
% % x_display_T(x_display_T<(max(x_display_T(:))*0.3))=0; %to reduce the artefact
% figure(exp_dFMT.figHdl_ExpSetting),hold on
% patch(isosurface(X,Y,Z,x_display_T,max(x_display_T(:))*iso_int), ...
%     'facecolor', '#00ffff', 'edgecolor', '#00ff00');
% hold off
