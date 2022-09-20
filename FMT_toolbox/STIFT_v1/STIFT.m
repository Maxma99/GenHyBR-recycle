%% information************************************************************
% Smart Toolkit for Fluorescence Tomography (STIFT)
% v 1.0
% Wuwei Ren, PhD
% 1 Institute of Biomedical Engineering, ETH Zurich
% 2 BORL, University Hospital Zurich
% Date: 
% Feb 5, 2019
%
%%************************************************************************
% 
% Discription:
% A new object-oriented version of FMT simulation/reconstruction with an 
% interactive GUI. It gives the possibility to use a class called STIFT. This
% class interacts with another two classes developed by others: TOAST++ and
% iso2Mesh, namely. Each FMT experiment can be customized with assist of STIFT
%
% General procedure:
% Step 1: Mesh creation
% Step 2: Optical coefficient and fluorescence assignment
% Step 3: Sources and detector setting
% Step 4: Running excitation process (excitation map recorded)
% Step 5: Running emission process (emission map recorded)
% Step 6: Reconstruction
% * the concept of Finite State Machine (FSM) is applied in this version
% * STIFT is working in 5 modules: 'configuration', 'simulator',
% * 'controller', 'reconstructor', 'viewer'.
% machine_mode: the state/mode of machine 
% 1 - Virtual FMT measurement&reconstruction
% 2 - Experimental FMT measurement&reconstruction
% 3 - Experimental FMT reconstruction
% *************************************************************************
% structure of the class STIFT
% module 1 'configuration': 
    % part 1.1 project info    
    % part 1.2 mesh
    % part 1.3 optical property
    % part 1.4 fluorescnce
    % part 1.5 illumination
    % part 1.6 detector
% module 2 'simulator': 
% module 3 'controller': 
%  under construction 
% module 4 'reconstructor':
% module 5 'viewer':
% *************************************************************************

classdef STIFT  < handle
    properties
    %% module 1 'configuration': 
        % part 1.1 project info
        project_directory % folder to save the project
        data_buffer_directory % folder to save the data buffer
        measurement % measurement data (in mode 3)
        machine_mode % machine mode
        date % experiment date
        time % experiment time
        object_ID % the object ID (mouse ID)
        user_name % the name of user
        % part 1.2 mesh
        mesh_method
        mesh_finess
        toast_mesh
        phantom_parameter
        n_node
        n_element
        node
        node_region
        element
        eltp
        face
        grd
        % part 1.3 optical property
        optiProp
        % part 1.4 fluorescnce
        fluo
        % part 1.5 illumination
        laser
        qvec_method
        qvec
        % part 1.6 detector
        detector
        detector_noise
        mvec_method
        mvec        
    %% module 2 'simulator': 
        phi_Ex % excitation fluence
        det_Ex % excitation map
        det_Ex_recon % excitation map for reconstruction
        phi_Em % emission fluence
        det_Em % emission map   
    %% module 3 'controller': 
    %  under construction 
    %% module 4 'reconstructor':
        recon_grd
        recon_coarseness
        mea_mask_option
        mea_mask
        mea_mask_thr
        sol_mask_option
        sol_mask
        sol_mask_grd
        dispose_mea_mask_option
        dispose_mea_mask
        dispose_mea_mask_thr
        weighting_Matrix
        measure_array
        reconSetting
        fluorecon
        system_K
    %% module 5 'viewer':
        figHdl_switch
        figHdl_ExpSetting
        figHdl_Ex
        figHdl_Em
        figHdl_Recon
        figHdl_FluoSetting
    end

    methods
    %% module 1 'configuration': 
        function configurationSTIFT(obj) % v1
            obj=configurationSTIFT_1(obj);
        end
        function confiSimuSTIFT(obj) % v1
            obj=confiSimuSTIFT_1(obj);
        end
        function confiExpSTIFT(obj) % v1
            obj=confiExpSTIFT_1(obj);
        end
        % part 1.1 project info
        function projectInfoSTIFT(obj) % v1
            obj=projectInfoSTIFT_1(obj);
        end        
        % part 1.2 mesh
        function geneMeshSTIFT(obj) % v1
            obj=geneMeshSTIFT_1(obj);
        end
        % part 1.3 optical property
        function optiAssignSTIFT(obj,Grid_Mesh_basis) % v1
            obj = optiAssignSTIFT_1(obj,Grid_Mesh_basis);
        end
        % part 1.4 fluorescnce
        function fluoSettingSTIFT(obj,Grid_Mesh_basis) % v1
            obj = fluoSettingSTIFT_1(obj,Grid_Mesh_basis);
        end
        % part 1.5 illumination
        function laserSettingSTIFT(obj) % v1
            obj = laserSettingSTIFT_1(obj);
        end
        % part 1.6 detector
        function detectorSettingSTIFT(obj) % v1
            obj = detectorSettingSTIFT_1(obj);
        end
        % part special
        function registSTIFT(obj) % v1
            obj = registSTIFT_1(obj);
        end
    %% module 2 'simulator': 
        function calcuQvecSTIFT(obj)  % v1
            obj = calcuQvecSTIFT_1(obj);
        end
        function calcuMvecSTIFT(obj)  % v1
            obj = calcuMvecSTIFT_1(obj);
        end
        function calcuSystemKSTIFT(obj) % v1
            obj = calcuSystemKSTIFT_1(obj);
        end
        function exRunSTIFT(obj,system_K_file) % v1
            obj = exRunSTIFT_1(obj,system_K_file);
        end
        function emRunSTIFT(obj,system_K_file) % v1
            obj = emRunSTIFT_1(obj,system_K_file);
        end
    %% module 3 'controller': 
    %  under construction
        function controllerSTIFT(obj) % v1
            obj=controllerSTIFT_1(obj);
        end
    %% module 4 'reconstructor':
        function reconSTIFT(obj) % v1
            obj = reconSTIFT_1(obj);
        end
       function geneMeaMaskSTIFT(obj) % v1
            obj = geneMeaMaskSTIFT_1(obj);
       end
        function geneDisposeMeaMaskSTIFT(obj) % v1
            obj = geneDisposeMeaMaskSTIFT_1(obj);
        end
        function geneSolMaskSTIFT(obj) % v1
            obj = geneSolMaskSTIFT_1(obj);
        end
        function geneWMSTIFT(obj,system_K_file) % v1
            obj = geneWMSTIFT_1(obj,system_K_file); 
        end
        function inversSTIFT(obj)
            obj = inversSTIFT_1(obj); % v1
        end
        function loadMeasurement(obj)
            obj = loadMeasurement_1(obj); % v1
        end
    %% module 5 'viewer':
        function drawExpSettingSTIFT(obj) % v1
            obj = drawExpSettingSTIFT_1(obj);
        end
        function drawExDetSTIFT(obj) % v1
            obj = drawExDetSTIFT_1(obj);
        end
        function drawEmDetSTIFT(obj) % v1
            obj = drawEmDetSTIFT_1(obj);
        end
        function drawReconSTIFT(obj) % v1
            obj = drawReconSTIFT_1(obj);
        end
        function drawFluoSettingSTIFT(obj) % v1
            obj = drawFluoSettingSTIFT_1(obj);
        end       
        function refineMeshSTIFT(obj) % v1
            obj = refineMeshSTIFT_1(obj);
        end
    end
end
