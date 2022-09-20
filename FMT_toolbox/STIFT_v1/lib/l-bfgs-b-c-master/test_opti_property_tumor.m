%% Pre-processing the MRI images: segmentation
% Wuwei Ren, 2017 May 30

close all;
clear all;
mouse_number=1;
level=0.1; % threshold for height map
xy_resolution = 20/128; % MRI xy resolution,mm
z_resolution = 0.5; % MRI z resolution, mm
x_slice_n = 128;
y_slice_n = 128;
z_slice_n = 35;
z_slice_n_update = round(z_slice_n*z_resolution/xy_resolution);

%% select the MRI raw images folder
fprintf(1,'loading MRI raw data ...\n');
Folder_name = uigetdir;
if isequal(Folder_name,0)
    disp('User selected Cancel');
    return;
end

%% read the raw images
raw_mri_data = zeros(x_slice_n,y_slice_n,z_slice_n);
for i= 1:z_slice_n
    if i<10
        image_name=['\mrimage_00' num2str(i) '.tif'];
    elseif i<100
        image_name=['\mrimage_0' num2str(i) '.tif'];
    else
        image_name=['\mrimage_' num2str(i) '.tif'];
    end
    temp_img_1=imread([Folder_name image_name]);
    raw_mri_data(:,:,i)=temp_img_1(:,:,1);
end

%% interplation 
[Xq,Yq,Zq] = meshgrid(1:x_slice_n,1:y_slice_n,linspace(1,z_slice_n,z_slice_n_update));
int_mri_data = interp3(raw_mri_data,Xq,Yq,Zq);
int_mri_data = abs(round(int_mri_data));

%% segmentation
seg_mri_data = zeros(x_slice_n,y_slice_n,z_slice_n_update);
for i =1:z_slice_n_update
    filter_temp_img=wiener2(int_mri_data(:,:,i),[8 8]);
    filter_temp_img_uint = uint8(filter_temp_img);
    BW = im2bw(filter_temp_img_uint, level);
    for j=1:y_slice_n
        line_BW = BW(:,j);
        if sum(line_BW)>0
           index_select = find(line_BW);
           range_select = min(index_select):max(index_select); 
           %range_select = min(index_select):110; 
           seg_mri_data(range_select,j,i) = 1;
        end
    end
    tmp_img_seg = seg_mri_data(:,:,i);
    edge_BW = edge(tmp_img_seg);
    se = strel('disk',1);
    edge_BW = imdilate(edge_BW,se);
    bone_BW = (filter_temp_img_uint>1)&(filter_temp_img_uint<25)&(tmp_img_seg>0);
    tmp_img_seg(bone_BW) = 3;
    tmp_img_seg(edge_BW) = 2;
    seg_mri_data(:,:,i) = tmp_img_seg;
    %{
    figure,
    subplot(1,3,1),imagesc(int_mri_data(:,:,i));
    subplot(1,3,2),imagesc(filter_temp_img);
    subplot(1,3,3),imagesc(BW);
    %}
end
int_mri_data_crop = int_mri_data.*(seg_mri_data>0);


%% 256*256*256 format saving for segmentation
final_mri_256 = zeros(256,256,256);
x_range_256 = 128-round(x_slice_n/2):(127-round(x_slice_n/2)+x_slice_n);
y_range_256 = 128-round(y_slice_n/2):(127-round(y_slice_n/2)+y_slice_n);
z_range_256 = 128-round(z_slice_n_update/2):(127-round(z_slice_n_update/2)+z_slice_n_update);
final_mri_256(x_range_256,y_range_256,z_range_256) = seg_mri_data;

%% 256*256*256 format saving for anatomy
final_mri_anatomy = zeros(256,256,256);
final_mri_anatomy(x_range_256,y_range_256,z_range_256) = int_mri_data_crop;

%% show the figures
%figure,imshow3Dfull(raw_mri_data);
%figure,imshow3Dfull(int_mri_data);
figure,imshow3Dfull(final_mri_256);
figure,imshow3Dfull(final_mri_anatomy);
%colormap(hot);
%figure,imshow3Dfull(int_mri_data_crop);
save([Folder_name '\tumor_seg.mat'],'final_mri_256');
save([Folder_name '\tumor_anatomy.mat'],'final_mri_anatomy');

