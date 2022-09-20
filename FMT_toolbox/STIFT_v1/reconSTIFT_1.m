% v1 
function obj = reconSTIFT_1(obj)
%% function illustration
% step 1: masking 
% step 2: generating a weighting matrix
% step 3: inversion

%% step 1: masking
    obj.geneMeaMaskSTIFT;
    obj.geneSolMaskSTIFT;
%% step 2: generating a weighting matrix
    obj.geneWMSTIFT(obj.system_K);
%% step 3: inversion
    obj.inversSTIFT;
end