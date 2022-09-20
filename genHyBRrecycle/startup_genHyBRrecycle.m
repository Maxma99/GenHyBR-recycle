%% startup_genHyBRrecycle.m
%
% Startup file for the HyBR recycling codes, AIRTools and IRtools needed.

directory = pwd;
path(directory, path)

path([directory, '/AIRToolsII-master'], path)
AIRToolsII_setup

path([directory, '/IRtools-master'], path)
IRtools_setup

path([directory, '/genHyBR'], path)
path([directory, '/HyBR'], path)
path([directory, '/genHyBRrecycle'], path)
path([directory, '/toeplitz'], path)
path([directory, '/seismic'], path)
clear directory




