%% v 1
%% module 1 'configuration': 
    % part 1.1 project info    
    % part 1.2 mesh
    % part 1.3 optical property
    % part 1.4 fluorescnce
    % part 1.5 illumination
    % part 1.6 detector

%%
function  obj  = configurationSTIFT_1(obj)
    switch obj.machine_mode
        case 1
            obj.confiSimuSTIFT;
        case 2
        case 3
            obj.confiExpSTIFT;
    end

end