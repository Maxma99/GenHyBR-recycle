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
close all;
%%
exp_dFMT.inversSTIFT;
exp_dFMT.drawExpSettingSTIFT;
exp_dFMT.drawExDetSTIFT;
exp_dFMT.drawEmDetSTIFT;
exp_dFMT.drawFluoSettingSTIFT;
exp_dFMT.drawReconSTIFT;

