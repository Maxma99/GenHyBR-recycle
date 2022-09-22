%% Instruction:
% a new STIFT protocol

%% clear memory
clearvars;
close all; 

%% Configurator
exp_dFMT = STIFT;
exp_dFMT.machine_mode = 1;% machine mode
exp_dFMT.configurationSTIFT;
% drawFroPhantom_new(exp_FMT.phantom_parameter.dim,exp_FMT.phantom_parameter.dl);

%% Simulator
exp_dFMT.calcuSystemKSTIFT;
exp_dFMT.exRunSTIFT(exp_dFMT.system_K);
exp_dFMT.emRunSTIFT(exp_dFMT.system_K); 

%% System Builder
exp_dFMT.geneMeaMaskSTIFT;
exp_dFMT.geneSolMaskSTIFT;
exp_dFMT.geneWMSTIFT(exp_dFMT.system_K);

%% Solver
exp_dFMT.inversSTIFT;

%% Viewer
exp_dFMT.drawExpSettingSTIFT;
exp_dFMT.drawExDetSTIFT;
exp_dFMT.drawEmDetSTIFT;
exp_dFMT.drawFluoSettingSTIFT;                                                                                                                                                                                                                                       
exp_dFMT.drawReconSTIFT;

