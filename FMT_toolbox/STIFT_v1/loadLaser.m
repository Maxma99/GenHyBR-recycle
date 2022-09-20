% v1
% p_list: giving the index numbers based on mat, real scale needed
% n_laser: number of laser points detected
function [p_list,n_laser] = loadLaser(thre,default_laser_folder,laser_used)
    tic
    fprintf(1,'loading laser from registered excitation map... ');
    [laser_ex_file,laser_ex_folder] = uigetfile('*.mat','laser from registered excitation map ',...
        default_laser_folder);
    if isequal(laser_ex_file,0)
        disp('no excitation selected')
        return
    else
        fprintf(1,'(found) \n');
    end
    laser_index = zeros(length(laser_used),2);
    vari_ex = load([laser_ex_folder laser_ex_file]);
    fn = fieldnames(vari_ex);
    mat_ex = vari_ex.(fn{1});    
    dim_mat_ex = size(mat_ex);
    background_img = mat_ex(:,:,1);
    sum_img = zeros(dim_mat_ex(1),dim_mat_ex(2));
    % read and analyze the images from the stack
    for i = 1:length(laser_used)
        tmp_img = mat_ex(:,:,(1+laser_used(i)));
        tmp_img = double(tmp_img-background_img);
        %figure,imagesc(tmp_img);
        sum_img = sum_img + tmp_img;
        tmp_img_bw = (tmp_img>thre*max(tmp_img(:)));
        se = strel('disk',2);%default = 5
        tmp_img_bw = imclose(tmp_img_bw,se);
        tmp_img_bw = imopen(tmp_img_bw,se);
        %figure,imagesc(tmp_img_bw);
        if sum(tmp_img_bw(:)) >0
            stat_img = regionprops(tmp_img_bw,'centroid'); %每次一个点光源，得到一幅excitation图片.每次求重心坐标
            temp_array = stat_img.Centroid;
            if isempty(temp_array)
               laser_index(i,:) = 9999;
            else 
               laser_index(i,:) = stat_img.Centroid;
            end
        else
            laser_index(i,:) = 9999;
        end
    end
    
    figure,imagesc(sum_img);
    axis image
    title('Illumination points (detected)')
    colormap(gray);
    hold on
    for i = 1:length(laser_used)
        plot(laser_index(i,1),laser_index(i,2),'ro'); %laser_index：25*2 每个点的坐标
    end
    hold off
%     %%%%%%%%%%%%%%%%%%% phantom 1 transmission
%     laser_index(5,:) = [];
%     laser_index(20,:) = [];
%     laser_index(23,:) = [];
    %%%%%%%%%%%%%%%%%%%
    p_list =laser_index;
    n_laser = size(p_list,1);
    fprintf(1,'done\n');
    toc
end