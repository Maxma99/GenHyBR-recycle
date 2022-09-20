% v1
% optical coeicient setting for two layers
function obj = optiAssignSTIFT_1(obj,Grid_Mesh_basis)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %         mua [real column vector]:
    %             nodal absorption coefficient [1/mm]
    %         mus [real column vector]:
    %             nodal reduced scattering coefficient [1/mm]
    %         ref [real column vector]:
    %             nodal refractive index
    %         freq [real]:
    %             modulation frequency [MHz]
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %       calculation of optical properties
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % reduced scattering coefficients: 1/mm
    %w_l = obj.laser.l_wavelength;% laser wavelength(nm)
    w_l = 670;% laser wavelength(nm)
    w_l_500 = w_l/500; % laser wavelength(nm) / 500nm
    
    a = [4.6 2.42 2.29 1.89]; 
    b = [1.421 1.611 0.716 1.286];
    mus_skin = a(1)*(w_l_500^(-b(1)));
    mus_brain = a(2)*(w_l_500^(-b(2)));
    mus_bone = a(3)*(w_l_500^(-b(3)));
    mus_softtissue = a(4)*(w_l_500^(-b(4)));
    
    % absorption: 1/mm
    mua_skin = 0.02;
    mua_brain = 0.015;
    mua_bone = 0.005;
    mua_softtissue = 0.02;    
    
    % refractive index: no unit
    ref_skin = 1.4;
    ref_brain = 1.4;
    ref_bone = 1.4;
    ref_softtissue = 1.4;
    
    % phantom properties
    % number 5
%     mua_A = 0.007;% mua in material A
%     mus_A = 0.87;% mus in material A
%     ref_A = 1.44;% ref in material A
    mua_A = 0.0052;% mua in material A
    mus_A = 0.99;% mus in material A
    ref_A = 1.44;% ref in material A
    
    % phantom properties
    % number 6
    mua_B = 0.0312;% mua in material B
    mus_B = 0.99;% mus in material B
    ref_B = 1.44;% ref in material B
    
    % phantom properties
    % number 7
    mua_C = 0.02;% mua in material A
    mus_C = 1;% mus in material A
    ref_C = 1.4;% ref in material A   
    
    % phantom properties
    % number 8
    mua_D = 0.0053;% mua in material A
    mus_D = 0.97 ;% mus in material A
    ref_D = 1.4;% ref in material A   
    
    % phantom properties
    % number 8
    mua_E =  0.033;% mua in material A
    mus_E = 0.89 ;% mus in material A
    ref_E = 1.4;% ref in material A 
    
    switch obj.optiProp.method
        case 'regu' 
            %digimouse
%             o_mua = zeros(obj.n_node,1);
%             o_mus = zeros(obj.n_node,1);
%             o_ref = zeros(obj.n_node,1);
%             for j = 1:obj.n_node
%                 if(obj.node_region(j)~=0)
%                     o_mua(j) = obj.optiProp.list_mua(obj.node_region(j));
%                     o_mus(j) = obj.optiProp.list_mus(obj.node_region(j));
%                     o_ref(j) = obj.optiProp.list_ref(obj.node_region(j));
%                 end
%             end
%             obj.optiProp.o_mua  = o_mua;
%             obj.optiProp.o_mus  = o_mus;
%             obj.optiProp.o_ref  = o_ref;
            % simu
            mua_g = drawBoxPhantom_new(obj.phantom_parameter.dim, ...
                obj.phantom_parameter.dl, ...
                obj.phantom_parameter.i_pos, ...
                obj.phantom_parameter.i_dim, ...
                obj.optiProp.list_mua);
            mus_g = drawBoxPhantom_new(obj.phantom_parameter.dim, ...
                obj.phantom_parameter.dl, ...
                obj.phantom_parameter.i_pos, ...
                obj.phantom_parameter.i_dim, ...
                obj.optiProp.list_mus);
            ref_g = drawBoxPhantom_new(obj.phantom_parameter.dim, ...
                obj.phantom_parameter.dl, ...
                obj.phantom_parameter.i_pos, ...
                obj.phantom_parameter.i_dim, ...
                obj.optiProp.list_ref);
            obj.optiProp.o_mua  = Grid_Mesh_basis.Map ('B->M', reshape(mua_g,[],1));
            %figure,imshow3Dfull(mua_g);
            obj.optiProp.o_mus  = Grid_Mesh_basis.Map ('B->M', reshape(mus_g,[],1));
            obj.optiProp.o_ref  = Grid_Mesh_basis.Map ('B->M', reshape(ref_g,[],1));            
        case 'load'
            temp_load = load([obj.optiProp.folder obj.optiProp.file]);
            fn = fieldnames(temp_load);
            croped_segment_img = temp_load.(fn{1});

            %% mua
            mua_g = croped_segment_img;
            mua_g(croped_segment_img == 1) = mua_softtissue;
            mua_g(croped_segment_img == 2) = mua_skin;
            mua_g(croped_segment_img == 3) = mua_bone;
            mua_g(croped_segment_img == 4) = mua_brain;
            mua_g(croped_segment_img == 5) = mua_A;
            mua_g(croped_segment_img == 6) = mua_B;
            mua_g(croped_segment_img == 7) = mua_C;
            mua_g(croped_segment_img == 8) = mua_D;
            mua_g(croped_segment_img == 9) = mua_E;
            mua_g(croped_segment_img == 0) = 0;
            %mua_g(croped_segment_img == 2) = mua_softtissue;
            %mua_g(croped_segment_img == 3) = mua_softtissue;
            %mua_g(croped_segment_img == 0) = mua_softtissue;

            %% mus
            mus_g = croped_segment_img;
            mus_g(croped_segment_img == 1) = mus_softtissue;
            mus_g(croped_segment_img == 2) = mus_skin;
            mus_g(croped_segment_img == 3) = mus_bone;
            mus_g(croped_segment_img == 4) = mus_brain;
            mus_g(croped_segment_img == 5) = mus_A;
            mus_g(croped_segment_img == 6) = mus_B;
            mus_g(croped_segment_img == 7) = mus_C;
            mus_g(croped_segment_img == 8) = mus_D;
            mus_g(croped_segment_img == 9) = mus_E;
            mus_g(croped_segment_img == 0) = 0;
            %mus_g(croped_segment_img == 2) = mus_softtissue;
            %mus_g(croped_segment_img == 3) = mus_softtissue;
            %mus_g(croped_segment_img == 0) = mus_softtissue;

            %% ref
            ref_g = croped_segment_img;
            ref_g(croped_segment_img == 1) = ref_softtissue;
            ref_g(croped_segment_img == 2) = ref_skin;
            ref_g(croped_segment_img == 3) = ref_bone;    
            ref_g(croped_segment_img == 4) = ref_brain;  
            ref_g(croped_segment_img == 5) = ref_A; 
            ref_g(croped_segment_img == 6) = ref_B; 
            ref_g(croped_segment_img == 7) = ref_C; 
            ref_g(croped_segment_img == 8) = ref_D;
            ref_g(croped_segment_img == 9) = ref_E;
            ref_g(croped_segment_img == 0) = 1.4;

            obj.optiProp.o_mua  = Grid_Mesh_basis.Map ('B->M', reshape(mua_g,[],1)); % B->M map to solution basis
            obj.optiProp.o_mua(obj.optiProp.o_mua<=0.001) = mua_skin;

            obj.optiProp.o_mus  = Grid_Mesh_basis.Map ('B->M', reshape(mus_g,[],1));
            obj.optiProp.o_mus(obj.optiProp.o_mus<=0.5) = mus_skin;

            obj.optiProp.o_ref  = Grid_Mesh_basis.Map ('B->M', reshape(ref_g,[],1));
            obj.optiProp.o_ref(obj.optiProp.o_ref<=1) = ref_skin;

            %figure,imshow3Dfull(mua_g);
            %colormap(jet)

            %figure,imshow3Dfull(mus_g);
            %colormap(jet)
            
    end