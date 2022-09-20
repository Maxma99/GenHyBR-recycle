%% v 1 configuration for experimental FMT
%% module 1 'configuration': 
    % part 1.1 project info    
    % part 1.2 mesh
    % part 1.3 optical property
    % part 1.4 fluorescnce (no need)
    % part 1.5 illumination
    % part 1.6 detector

%%
function  obj  = confiExpSTIFT_1(obj)
    % part 1.1 project info
    obj.projectInfoSTIFT;
    % part 1.2 mesh + registration
    obj.mesh_method = 3;
    %    mesh_method: phantom selection
    %                 1 - standard geometry generated from iso2mesh
    %                 2 - standard geometry loaded from gmsh
    %                 3 - complex geometry generated from iso2mesh with
    %                 datasets
    %                 4 - mouse atlas
    obj.phantom_parameter.setup = 1; % a phantom setup switch 1/0
    if obj.phantom_parameter.setup == 1 % build up a slab phantom
        obj.phantom_parameter.dx = 1; % mm
        obj.phantom_parameter.dy = obj.phantom_parameter.dx; % mm
        obj.phantom_parameter.dz = obj.phantom_parameter.dx; % mm
        obj.phantom_parameter.x = 54; % mm 15
        obj.phantom_parameter.y = 54; % mm 30 
        obj.phantom_parameter.z = 14; % mm 60
        obj.phantom_parameter.cap_d = 3; % mm
        obj.phantom_parameter.cap_size = 1; %mm
        obj.phantom_parameter.dim = [obj.phantom_parameter.x obj.phantom_parameter.y obj.phantom_parameter.z]';
        obj.phantom_parameter.dl=[obj.phantom_parameter.dx obj.phantom_parameter.dy obj.phantom_parameter.dz]';
        obj.phantom_parameter.i_pos = [obj.phantom_parameter.dim/2 ...
            obj.phantom_parameter.dim/4]; % inclusion positions
        obj.phantom_parameter.i_dim = [obj.phantom_parameter.dim/4 ...
            obj.phantom_parameter.dim/8]; % inclusion dimensions
        obj.grd = ([obj.phantom_parameter.x obj.phantom_parameter.y obj.phantom_parameter.z]./...
            [obj.phantom_parameter.dx obj.phantom_parameter.dy obj.phantom_parameter.dz]) + 1;
    end
    obj.mesh_finess = 1;
    obj.geneMeshSTIFT;


% toast_mesh = toastMesh(newnode,newelem,eltp);
    obj.toast_mesh = toastMesh(obj.node,obj.element,obj.eltp);
    %Grid_Mesh_basis = toastBasis(obj.toast_mesh, obj.grd, 'CUBIC');
    obj.registSTIFT;
    
    % part 1.3 optical property
    %         optiProp.o_mua [real column vector]:
    %             nodal absorption coefficient [1/mm]
    %         optiProp.o_mus [real column vector]:
    %             nodal reduced scattering coefficient [1/mm]
    %         optiProp.o_ref [real column vector]:
    %             nodal refractive index
    %         optiProp.method [string]:
    %             'regu' or 'load'
    obj.optiProp.method = 'load'; % 'regu' or 'load'
    switch obj.optiProp.method
        case 'regu' % build the optical properties
            obj.optiProp.list_mua = [0.01 0.01 0.01]';
            obj.optiProp.list_mus = [1 1 1]';
            obj.optiProp.list_ref = [1.4 1.4 1.4]';
            obj.optiProp.opti_grd = obj.grd;
            opti_Grid_Mesh_basis = toastBasis(obj.toast_mesh,...
                obj.optiProp.opti_grd);
            obj.optiAssignSTIFT(opti_Grid_Mesh_basis);
        case 'load' % load the optical properties
            [obj.optiProp.file,obj.optiProp.folder] = uigetfile...
                ('*.mat','load optical property (e.g. cropped)',obj.project_directory);
            temp_load = load([obj.optiProp.folder obj.optiProp.file]);
            fn = fieldnames(temp_load);
            croped_segment_img = temp_load.(fn{1});
            obj.optiProp.opti_grd =size(croped_segment_img);
            clear croped_segment_img
            opti_Grid_Mesh_basis = toastBasis(obj.toast_mesh,...
                obj.optiProp.opti_grd); % solution basis
            obj.optiAssignSTIFT(opti_Grid_Mesh_basis);
    end

    % part 1.4 fluorescnce: skipped
    %
    % part 1.5 illumination
    obj.laser.l_design ='regu'; %'regu' or 'load'
    obj.qvec_method = 'gaussian-neumann';
    % 'regu': define the regular pattern
    % 'load': input from excitation map
    switch obj.laser.l_design    
        case 'regu'
            obj.laser.l_power = 1; % W/mm2
            obj.laser.l_width = 1; % mm
            obj.laser.l_freq = 0; % CW mode
            obj.laser.l_wavelength = 630; % nm
            obj.laser.l_scanGapTime = 1; % s
%             obj.laser.l_direction = [0 0 -1]; % laser direction
            obj.laser.l_direction = [0 0 1]; % laser direction
            obj.laser.l_edge = 21; % mm
            obj.laser.l_number = [10 10]; 
            obj.laser.l_pos_used = 1:prod(obj.laser.l_number); %number of points used
            if obj.laser.l_direction(3) ~= 0
                obj.laser.l_region= [min(obj.node(:,1))+obj.laser.l_edge ...
                    max(obj.node(:,1))-obj.laser.l_edge; ...
                    min(obj.node(:,2))+obj.laser.l_edge ...
                    max(obj.node(:,2))-obj.laser.l_edge]; 
                % 2D plane two corners, [x1 x2;y1 y2]   
            end
            if obj.laser.l_direction(2) ~= 0
            	obj.laser.l_region= [min(obj.node(:,1))+obj.laser.l_edge ...
                    max(obj.node(:,1))-obj.laser.l_edge; ...
                    min(obj.node(:,3))+obj.laser.l_edge ...
                    max(obj.node(:,3))-obj.laser.l_edge]; 
                % 2D plane two corners, [x1 x2;z1 z2]
            end                    
            if obj.laser.l_direction(1) ~= 0
            	obj.laser.l_region= [min(obj.node(:,2))+obj.laser.l_edge ...
                    max(obj.node(:,2))-obj.laser.l_edge; ...
                    min(obj.node(:,3))+obj.laser.l_edge ...
                    max(obj.node(:,3))-obj.laser.l_edge]; 
                % 2D plane two corners, [y1 y2;z1 z2]
            end
            obj.laserSettingSTIFT; % create the laser positio
        case 'load'
            obj.laser.l_power = 1.5; % W/mm2
            obj.laser.l_width = 1; % mm
            obj.laser.l_freq = 0; % CW mode
            obj.laser.l_wavelength = 630; % nm
            obj.laser.l_scanGapTime = 1; % s
            obj.laser.l_direction = [0 0 1]; % laser direction
            obj.laser.l_pos_used = 1:prod(obj.laser.l_number); %number of points used
            obj.laser.l_number = [10 10];% reflection ([7 7]) 
            obj.laserSettingSTIFT; % create the laser position
    end
        
    % part 1.6 detector
    obj.detector.d_number = [obj.phantom_parameter.x+1 obj.phantom_parameter.y+1];
    obj.detector.d_dim = [1 1]; % mm
    obj.detector.d_efficiency = 1;
    obj.detector.d_virtualPos = [];
    obj.detector.d_normal_dA = [0 0 -1];% normal dA direction, right hand rule
    obj.detector.d_edge = 0; % mm
    obj.detector.d_mode = 2; % detector mode:1 & 2, only useful for MRI input
    obj.detector.d_region = [];
    obj.mvec_method = 'gaussian';
    obj.detector_noise.option = 0; % noise switch 0 or 1
    obj.detector_noise.level = 0; % noise level; normally nonnegative value; level<5
    obj.detectorSettingSTIFT; % create the detector position
    % load both detector and laser
    fprintf(1,'Laser and Detector are loaded and transfered to TOAST...\n');
    tic
    obj.toast_mesh.SetQM(obj.laser.l_position,obj.detector.d_virtualPos); % SetQM 方法将源位置和探测器位置附加到网格
    fprintf(1,'done\n');
    toc
    % plot网格以及source和detector位置
%     obj.toast_mesh.Display
%     hold on
%     plot3(obj.laser.l_position(:,1),obj.laser.l_position(:,2),obj.laser.l_position(:,3),'ro','MarkerFaceColor','r');
%     plot3(obj.detector.d_virtualPos(:,1),obj.detector.d_virtualPos(:,2),obj.detector.d_virtualPos(:,3),'bx','MarkerFaceColor','b');
%     hold off

    %% Reconstruction configuration
    obj.recon_coarseness = 2;
    obj.recon_grd = round((obj.grd -1)./obj.recon_coarseness) + 1;
    obj.mea_mask_option = 1;
    obj.mea_mask_thr = 2e-2;
    obj.mea_mask = [];
    obj.dispose_mea_mask_option = 1;
    obj.dispose_mea_mask_thr = 2e-2;
    obj.dispose_mea_mask = [];
    obj.sol_mask_option = 1;
    obj.sol_mask = [];
    obj.sol_mask_grd = [];
    obj.weighting_Matrix = [];
    obj.measure_array = [];
    obj.reconSetting.Method = 10;
    % reconstruction method options: 
    %   1,2: ART w/o & w/ regularization
    %   3,4: MATLAB simple backslash w/o & w/ regularization
    %   5,6: CG method: pcg w/o & w/ regularization
    %   7,8: lsqr w/o & w/ regularization
    %   9,10: l1 regularization using ISTA/FISTA
    %   11: l1 joint l2 
    %   12: l1 joint TV
    %   13, 14: non negative least square w/o & w/ regularization
    obj.reconSetting.iterTime = 500;
    obj.reconSetting.lambdaFactor = 1e-4;
    obj.reconSetting.tolerror = 1e-8;
    obj.reconSetting.recon_use_real_ex = 0;
    obj.fluorecon = [];
    %}
end