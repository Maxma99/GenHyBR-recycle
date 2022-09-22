% v1
function obj = emRunSTIFT_1(obj,system_K_file)    
    tic
    fprintf(1,'Emission simulation...\n');
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
   
    switch obj.optiProp.method
        case 'regu' % build the optical properties
            Q_fluo=obj.fluo.fluoDis_m*ones(1,length(obj.laser.l_position(:,1)));
            Em_qvec = Q_fluo.*obj.phi_Ex;
            obj.phi_Em = system_K\Em_qvec;
            obj.det_Em = mvec.' * obj.phi_Em;
            if obj.detector_noise.option
                det_size = size(obj.det_Em);
                det_noise = rand(det_size) .* mean(obj.det_Em(:)) .* obj.detector_noise.level;
                obj.det_Em = obj.det_Em + det_noise;
            end
            fprintf(1,'done\n');
        case 'load' % load the optical properties
            [em_file,em_folder] = uigetfile('*.mat','registered emission map ',...
                 obj.project_directory);
             map_em = load([em_folder em_file]);
             fn = fieldnames(map_em);
             mat_em = map_em.(fn{1}); 
%              mat_em = mat_em(:,:,2:end-1);
%              background_img = mat_em(:,:,1);
%              %%%%%%%%%%%%%%%%%%%%%%% phantom1 transmission
%              mat_em(:,:,5) = [];
%              %%%%%%%%%%%%%%%%%%%%%%%
%              obj.det_Em = zeros(obj.grd(1)*obj.grd(2),size(qvec,2));
%              mat_em_x_ds = ceil(resample(1:size(mat_em,1),obj.grd(1),size(mat_em,1)));
%              mat_em_y_ds = ceil(resample(1:size(mat_em,2),obj.grd(2),size(mat_em,2)));
%              mat_em_x_ds(54:61) = [339,346,353,359,366,372,379,386];
%              mat_em_y_ds(24:31) = [143,149,156,163,170,176,183,189];
%              for i = 1:size(mat_em,3)
%                  temp = mat_em(mat_em_x_ds,mat_em_y_ds,i)';
%                  obj.det_Em(:,i) = temp(:);
%              end
            for i = 1:size(mat_em,3)
                
%                 corr_em = corr_em - background_img;
%                 corr_em(corr_em<0) = 0;
                mat_em(:,:,i) =  mat_em(:,:,i)-min(mat_em(:,:,i));
                corr_em =  mat_em(:,:,i);
               ratex = size(mat_em,1)/max(obj.detector.d_virtualPos(:,1));
               ratey = size(mat_em,2)/max(obj.detector.d_virtualPos(:,2));
               a=ceil(ratex * obj.detector.d_virtualPos(:,1)); 
               a(a==0) = 1;
                b=ceil(ratey * obj.detector.d_virtualPos(:,2));
                 b(b==0) = 1;
                 %% circle
%                 temp = [];
%                  for j = 1:length(a)
% %                      tmp = corr_em(b(j):b(j),a(j));
%                      temp(j) = corr_em(b(j),a(j));
%                  end
%                  obj.det_Em(:,i) = temp;
                %% square
%                  temp = zeros(obj.grd(1)*obj.grd(2),1);
%                  for j = 1:length(a)
%                      if b(j)>ratex/2 && b(j)<size(corr_em,2)-ratex/2 &&...
%                              a(j)>ratey/2 && a(j)<size(corr_em,1)-ratey/2
%                          temp_corr_em = corr_em(a(j)-floor(ratex/2):a(j)+floor(ratex/2),b(j)-floor(ratex/2):b(j)+floor(ratex/2));
%                          temp(j) = mean(temp_corr_em(:));
%                      else
%                          temp(j) = corr_em(a(j),b(j));
%                      end
%                         
%                  end
%                  obj.det_Em(:,i) = temp;
                 a = reshape(a,obj.grd(2),obj.grd(1));
                 b = reshape(b,obj.grd(2),obj.grd(1));       
               temp = mat_em(a(1,:)',b(:,1),i)';
               obj.det_Em(:,i) = temp(:);
            end
            det_Em = obj.det_Em;
                                %%%%%%%%%%%%%%%%%%% phantom 1 transmission
%     obj.det_Em(:,5) = [];
%     obj.det_Em(:,20) = [];
%     obj.det_Em(:,23) = [];
            save([obj.data_buffer_directory '/det_Em.mat'],'det_Em');
    end
    toc
end