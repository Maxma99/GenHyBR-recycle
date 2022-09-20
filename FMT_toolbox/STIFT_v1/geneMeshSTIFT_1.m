function  obj  = geneMeshSTIFT_1(obj)
%    mesh_method: phantom selection
%                 1 - standard geometry generated from iso2mesh
%                 2 - standard geometry loaded from gmsh
%                 3 - complex geometry generated from iso2mesh with
%                 datasets
%                 4 - mouse atlas
    fprintf(1,'Generating a mesh...\n');
    tic
    switch obj.mesh_method
        case 1
            fprintf(1,'1 - standard geometry generated from iso2mesh\n');
            iso2mesh_method = 'cubic_2';
            switch iso2mesh_method
                case 'cubic_1'
                    start_point = [0 0 0];
                    end_point = obj.phantom_parameter.dim';
                    [obj.node,obj.element]=meshgrid5...
                        (start_point(1):obj.phantom_parameter.dx:end_point(1),...
                        start_point(2):obj.phantom_parameter.dy:end_point(2),...
                        start_point(3):obj.phantom_parameter.dz:end_point(3));
                    obj.face=volface(obj.element);                      
                case 'cubic_2'
                    start_point = [0 0 0];
                    end_point = obj.phantom_parameter.dim';
                    %meshabox函数在window和macos上运算结果不同
                    [obj.node,obj.face,obj.element] = ...
                        meshabox(start_point,end_point,obj.mesh_finess*(obj.phantom_parameter.dx)^3,1);
                case 'cubic_3'
                    start_point = [0 0 0];
                    end_point = obj.phantom_parameter.dim';
                    [obj.node,obj.element]=meshgrid6...
                        (start_point(1):obj.phantom_parameter.dx:end_point(1),...
                        start_point(2):obj.phantom_parameter.dy:end_point(2),...
                        start_point(3):obj.phantom_parameter.dz:end_point(3));
                    obj.face=volface(obj.element);                      
            end
            obj.eltp = ones(length(obj.element),1)*3;
            obj.n_node = length(obj.node(:,1));
            obj.n_element = length(obj.element(:,1));  
        case 2
            fprintf(1,'2 - standard geometry loaded from gmsh\n');
            [FileName,PathName] = uigetfile('*.msh','Select gmsh file');
            if isequal(FileName,0)
                disp('User selected Cancel');
            else
                N_n_ID  = 8; % standard value is 8
                N_n     = dlmread([PathName FileName],'',[N_n_ID 0 N_n_ID 0]);
                %node_id = dlmread([PathName FileName],'',[N_n_ID+1 0 N_n_ID+N_n 0]);
                nodes   = dlmread([PathName FileName],'',[N_n_ID+1 1 N_n_ID+N_n 3]);
                N_e     = dlmread([PathName FileName],'',[N_n_ID+N_n+3 0 N_n_ID+N_n+3 0]);
                %element_id = dlmread([PathName FileName],'',[N_n_ID+N_n+4  0 N_n_ID+N_n+3+N_e 0]);
                elements   = dlmread([PathName FileName],'',[N_n_ID+N_n+4  5 N_n_ID+N_n+3+N_e 8]);
            end
            gmsh_factor = 10; % standard value is 10
            % obj.grd = round((max(nodes)-min(nodes))*gmsh_factor);
            obj.node = nodes * gmsh_factor;
            obj.node = obj.node - repmat(min(obj.node),N_n,1);
            obj.element = elements;
            obj.face = volface(obj.element);
            obj.n_node = N_n;
            obj.n_element = N_e;            
            obj.eltp = ones(obj.n_element,1)*3;
        case 3
%             obj_temp = readObj('100000newpoints3D_orig.obj');
%             v = obj_temp.v;
%             f = obj_temp.f.v;
% 
%             % img = surf2volz(v,f,0:0.1:20,0:0.1:20,0:0.1:11);
%             % figure;imshow3Dfull(img)
% 
%             start_point = [0,0,0];
%             end_point = [54,54,14];
%             grd = end_point+1;
%             max_volume = 0.5;
%             nodesize = 1;
% 
%             opt = 0.2;
%             [v1,f1] = remeshsurf(v,f,opt);
% 
%             % [v1,f1] = meshcheckrepair(v,f,'deep');
%             % [v1,f1] = meshresample(v,f,0.2);
%             % figure;plotmesh(v1,f1)
% 
% 
%             [obj.node,obj.element] = surf2mesh(v1,f1,start_point,end_point,1,max_volume);
%             obj.element = obj.element(:,1:4);
% %             figure;plotmesh(newnode,newelem)
%             % figure;plotmesh(v1,f1)
% 
%             obj.face = volface(obj.element);
%             obj.eltp = ones(length(obj.element),1)*3;
% 
%             obj.n_node = length(obj.node(:,1));
%             obj.n_element = length(obj.element(:,1)); 
            
            
            
            dataset_dx = obj.measurement.pixel_size; % the real size of one pixel
            dataset_dy = obj.measurement.pixel_size; % the real size of one pixel
            dataset_dz = obj.measurement.pixel_size; % the real size of one pixel  
            % 3 - 使用数据集从 iso2mesh 生成的复杂几何
            fprintf(1,'3 - complex geometry generated from iso2mesh with datasets\n');
            [FileName,PathName] = uigetfile('*.mat','Select meshing file (from lib)');
            temp_load = load ([PathName FileName]);
            fn = fieldnames(temp_load); %返回该结构体的所有字段名称组成的元胞数组
            pre_full_phantom = temp_load.(fn{1});
            clear temp_load fn;
            full_phantom=fillholes3d(logical(pre_full_phantom>0),10);
            full_phantom=double(full_phantom);
            %figure,imshow3Dfull(full_phantom);
            %title('meshed object')
            phantom_dim = size(full_phantom);
            opt.radbound = 0.03*min(phantom_dim); % max radius of the Delaunay sphere 0.03 is a good value
            opt.distbound = opt.radbound/3; % maximum deviation from the specified isosurfaces
            opt.side='lower'; % - 'upper': threshold at upper interface; 'lower': threshold at lower interface
            opt.keepratio=0.8; % setting compression rate for each levelset
            opt.angbound = 45; % defines the miminum angle of a surface triangle
            dofix = 1; % 1: perform mesh validation&repair, 0: skip repairing
            m_method = 'cgalsurf';% 'cgalsurf' or 'cgalpoly':
            isovalues = 1; % isosurface value
            resample_ratio = 0.3 ;% to remove the dense corner
            iter = 2; % smoothing iteration number
            alpha = 0.8; % scaler, smoothing parameter, v(k+1)=alpha*v(k)+(1-alpha)*mean(neighbors)
            keepratio = 0.8; % 
            maxvol = obj.mesh_finess*0.5*(opt.radbound^3); % inside maxvol of element
            % generate the surface element first
            [no,face]= vol2surf(full_phantom,1:phantom_dim(1),...
                  1:phantom_dim(2),1:phantom_dim(3),opt,dofix,m_method,isovalues);
            [n_no,n_face] = meshresample(no,face(:,1:3),resample_ratio);
            [n_no,n_face] = meshcheckrepair(n_no,n_face);% correction
            s_n_no = sms(n_no,n_face(:,1:3),iter,alpha,'lowpass'); 
            [s_n_no,n_face] = meshcheckrepair(s_n_no,n_face);% correction
            [m_node,m_elem,m_face]=s2m(s_n_no,n_face(:,1:3),keepratio,maxvol);
            obj.node = m_node-repmat(min(m_node),length(m_node(:,1)),1);
            %obj.grd = round(max(obj.node)-min(obj.node));
            obj.node(:,1) = obj.node(:,1)*dataset_dx;
            obj.node(:,2) = obj.node(:,2)*dataset_dy;
            obj.node(:,3) = obj.node(:,3)*dataset_dz;
            obj.element = m_elem(:,1:4);
            obj.face = m_face(:,1:3);
            obj.n_node = length(obj.node(:,1));
            obj.n_element = length(obj.element(:,1));            
            obj.eltp = ones(obj.n_element,1)*3;
            %% to save the croped segmentation images
            crop_start = round(min(m_node));
            crop_end = round(max(m_node));
            range_x = crop_start(1):crop_end(1);
            range_y = crop_start(2):crop_end(2);
            range_z = crop_start(3):crop_end(3);
            cropped_segment_img = pre_full_phantom(range_x,range_y,range_z);
            save([obj.project_directory '\cropped_segment_img.mat'], 'cropped_segment_img');
            %figure,imshow3Dfull(cropped_segment_img);
            %title('cropped segmented image')
            %clear croped_segment_img
        case 4
            fprintf(1,'4 - mouse atlas\n');
            load('D:\1 code\01_project_s_decom\STIFT_v1\digimouse_mesh\Digimouse_Mesh_1L.mat');
            obj.face = face(:,1:3);
            obj.node = node;
            obj.element = elem(:,1:4);
            obj.n_element = length(obj.element(:,1));
            obj.eltp = ones(obj.n_element,1)*3;
            obj.n_node = length(obj.node(:,1));
            obj.node_region = zeros(obj.n_node,1);
            for j = 1:max(elem(:,end))
                temp = elem(elem(:,end)==j);
                temp = unique(temp);
                obj.node_region(temp)=j;
            end
            
    end
    fprintf(1,'done\n');
    toc
    
    %% generate mesh topology
%     if obj.machine_mode == 3
%         tic 
%         fprintf(1,'generating mesh topology...\n');
%         switch obj.mesh_method
%             case {1,2}
%                 max_node = max(obj.node);
%                 mesh_top = ones(round(max_node(1)/obj.measurement.pixel_size),...
%                     round(max_node(2)/obj.measurement.pixel_size));
%                 %figure,imagesc(mesh_top);
%                 %title('mesh surface');
%                 save([obj.project_directory '\mesh_top.mat'], 'mesh_top');
%             case 3
%                 mesh_im_reg_full = sum(full_phantom,3);
%                 mesh_im_reg_full_s = shiftdim(mesh_im_reg_full);
%                 %figure, imagesc(mesh_im_reg_full_s);
%                 max_m_node = round(max(m_node));
%                 min_m_node = round(min(m_node));
%                 mesh_top = mesh_im_reg_full_s(min_m_node(1):max_m_node(1),...
%                     min_m_node(2):max_m_node(2));
%                 %figure,imagesc(mesh_top);
%                 %title('mesh surface');
%                 save([obj.project_directory '\mesh_top.mat'], 'mesh_top');
%             case 4
%         end
%         fprintf(1,'done\n');
%         toc       
%     end
end