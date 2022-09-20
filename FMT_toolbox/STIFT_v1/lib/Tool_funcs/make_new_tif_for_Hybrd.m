%% Pre-processing for A2 phantoms for hybrid FMT/MRI system: make new tiff stack
% Wuwei Ren, 2017 Jun 26
% category 1: homogeneous slab phanom 

clear all
close all
folder_name ='C:\Users\renw\Desktop\data\201705Hybrid_FMT_MRI_Elmer\4_FMT with Phantom reflection and transmission\';
wt_file_name = 'whitelight.tif';
ex_file_name = 'REFLECTIOn_excitation.tif';
em_file_name = 'REFLECTION_emission.tif';
cropped_window_x = 514:892;
cropped_window_y = 300:602;
x_l = length(cropped_window_x);
y_l = length(cropped_window_y);
raw_size = max(x_l,y_l);
wt_raw = imread([folder_name wt_file_name]);
wt_cropped = zeros(raw_size);
wt_cropped(1:x_l,1:y_l) = wt_raw(cropped_window_x,cropped_window_y);
wt_cropped_256 =  uint16(imresize(wt_cropped,[256 256]));
%% write white light
imwrite(wt_cropped_256,'wt_256.tif');
figure,imagesc(wt_raw);
axis equal
figure,imagesc(wt_cropped);
figure,imagesc(wt_cropped_256);



%% write excitation map
%temp_256 = zeros(256);
info = imfinfo([folder_name ex_file_name]);
num_images = numel(info);
for k = 1:num_images
    temp_cropped = zeros(raw_size);
    A = imread([folder_name ex_file_name], k);
    temp_cropped(1:x_l,1:y_l) = A(cropped_window_x,cropped_window_y);
    temp_256 = uint16(imresize(temp_cropped,[256 256]));
    
    if k ==1
        imwrite(temp_256, 'ex_stack.tif');
    else
        imwrite(temp_256, 'ex_stack.tif', 'writemode', 'append');
    end
    
end

%% write emission map
%temp_256 = zeros(256);
info = imfinfo([folder_name em_file_name]);
num_images = numel(info);
for k = 1:num_images
    temp_cropped = zeros(raw_size);
    A = imread([folder_name em_file_name], k);
    temp_cropped(1:x_l,1:y_l) = A(cropped_window_x,cropped_window_y);
    temp_256 = uint16(imresize(temp_cropped,[256 256]));
    
    if k ==1
        imwrite(temp_256, 'em_stack.tif');
    else
        imwrite(temp_256, 'em_stack.tif', 'writemode', 'append');
    end
    
end

