%% Instruction:
% a new STIFT protocol

%% clear memory
clearvars;
close all;
%% Configurator
exp_FMT = STIFT;
exp_FMT.machine_mode = 3;% machine mode
exp_FMT.configurationSTIFT;

%% Simulator
exp_FMT.calcuSystemKSTIFT;
exp_FMT.exRunSTIFT(exp_FMT.system_K);
%exp_FMT.emRunSTIFT(exp_FMT.system_K); % skipped

%% Reconstructor
exp_FMT.geneMeaMaskSTIFT;
exp_FMT.geneSolMaskSTIFT;
exp_FMT.geneWMSTIFT(exp_FMT.system_K);
exp_FMT.loadMeasurement;
exp_FMT.inversSTIFT;

%% Viewer
exp_FMT.drawExpSettingSTIFT;
exp_FMT.drawExDetSTIFT;
exp_FMT.drawEmDetSTIFT;
%exp_FMT.drawFluoSettingSTIFT;
exp_FMT.drawReconSTIFT;

%% save the reconstruction result
xx = exp_FMT.fluorecon;
save([exp_FMT.project_directory '\recon.mat'],'xx'); 