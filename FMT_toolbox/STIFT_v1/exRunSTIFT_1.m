% v1
function obj = exRunSTIFT_1(obj,system_K_file)
    tic
    fprintf(1,'Excitation simulation...\n');
    %% load system_K
    system_K_s = load(system_K_file);
    system_K = system_K_s.system_K;
    clear system_K_s
    %%%%%
    %% load mvec
    mvec_s = load(obj.mvec);
    mvec = mvec_s.mvec;
    clear mvec_s
    %% load qvec
    qvec_s = load(obj.qvec);
    qvec = qvec_s.qvec;
    clear mvec_s
    %%%%% 
    obj.phi_Ex = system_K\qvec;
    switch obj.optiProp.method
        case 'regu' % build the optical properties
            
            obj.det_Ex = mvec.' * obj.phi_Ex; % adding noise is also possible
            obj.det_Ex_recon = obj.det_Ex; % save for reconstruction
            if obj.detector_noise.option
                det_size = size(obj.det_Ex);
                det_noise = rand(det_size) .* mean(obj.det_Ex(:)) .* obj.detector_noise.level;
                obj.det_Ex = obj.det_Ex + det_noise;
            end
            fprintf(1,'done\n');
        case 'load' % load the optical properties
%             obj.det_Ex = mvec.' * obj.phi_Ex; % 表面上可测量到的data
 [ex_file,ex_folder] = uigetfile('*.mat','registered excitation map ',...
                 obj.project_directory);
             map_ex = load([ex_folder ex_file]);
             fn = fieldnames(map_ex);
             mat_ex = map_ex.(fn{1}); 

            for i = 1:size(mat_ex,3)
                corr_ex =  mat_ex(:,:,i);
                mat_ex(:,:,i) =  mat_ex(:,:,i)-min(mat_ex(:,:,i));
               ratex = size(mat_ex,1)/max(obj.detector.d_virtualPos(:,1));
               ratey = size(mat_ex,2)/max(obj.detector.d_virtualPos(:,2));
               a=ceil(ratex * obj.detector.d_virtualPos(:,1)); 
               a(a==0) = 1;
                b=ceil(ratey * obj.detector.d_virtualPos(:,2));
                 b(b==0) = 1;

                 a = reshape(a,obj.grd(2),obj.grd(1));
                 b = reshape(b,obj.grd(2),obj.grd(1));
               temp = mat_ex(a(1,:)',b(:,1),i)';
               obj.det_Ex(:,i) = temp(:);
               
            end
    end

    det_Ex = obj.det_Ex;
    save([obj.data_buffer_directory '/det_Ex.mat'],'det_Ex');
end