%% Instruction:
% a new STIFT protocol

%% clear memory
clearvars;
close all; 

%% Configurator
exp_FMT = STIFT;
exp_FMT.machine_mode = 1;% machine mode
exp_FMT.configurationSTIFT;
% drawFroPhantom_new(exp_FMT.phantom_parameter.dim,exp_FMT.phantom_parameter.dl);

%% Simulator
exp_FMT.calcuSystemKSTIFT;
exp_FMT.exRunSTIFT(exp_FMT.system_K);
exp_FMT.emRunSTIFT(exp_FMT.system_K); 

%% Reconstructor
exp_FMT.geneMeaMaskSTIFT;
exp_FMT.geneSolMaskSTIFT;
exp_FMT.geneWMSTIFT(exp_FMT.system_K);
exp_FMT.inversSTIFT;

%% Viewer
exp_FMT.drawExpSettingSTIFT;
exp_FMT.drawExDetSTIFT;
exp_FMT.drawEmDetSTIFT;
exp_FMT.drawFluoSettingSTIFT;                                                                                                                                                                                                                                       
exp_FMT.drawReconSTIFT;

%%

% x = 1:6:330;					% 原始矩阵插值后的x轴坐标位置,即第1、5、9个    
% y = 1:6:330;					% y轴坐标含义同上
% [X,Y] = meshgrid(x,y);		% 过渡
% 
% x = 1:330;					% 插值后矩阵的x轴数值，即[1 2 3 4 5 6 7 8 9]
% y = 1:330;					% y轴坐标含义同上
% [Xq,Yq] = meshgrid(x,y);	% 过渡
% 
% B = interp2(X,Y,gt,Xq,Yq); 	% 生成线性插值后的矩阵
% B = interp2(X,Y,A,Xq,Yq); 	% 生成线性插值后的矩阵

% exp_FMT.fluorecon = l1FISTA;
% exp_FMT.drawExpSettingSTIFT;
% exp_FMT.drawReconSTIFT;
% exp_FMT.fluorecon = l1l2;
% exp_FMT.drawExpSettingSTIFT;
% exp_FMT.drawReconSTIFT;
% exp_FMT.fluorecon = l2CG;
% exp_FMT.drawExpSettingSTIFT;
% exp_FMT.drawReconSTIFT;
% exp_FMT.fluorecon = s1;
% exp_FMT.drawExpSettingSTIFT;
% exp_FMT.drawReconSTIFT;
% exp_FMT.fluorecon = s2;
% exp_FMT.drawExpSettingSTIFT;
% exp_FMT.drawReconSTIFT;
% exp_FMT.fluorecon = s1_gt;
% exp_FMT.drawExpSettingSTIFT;
% exp_FMT.drawReconSTIFT;