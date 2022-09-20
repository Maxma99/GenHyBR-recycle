%% Pre-processing for A2 phantoms for hybrid FMT/MRI system
% Wuwei Ren, 2017 Jun 26
% category 1: homogeneous slab phanom 

close all;
clear all;
category_phantom = 1;
xyz_resolution = 20/128; % MRI xy resolution,mm
x_slice_n = 256;
y_slice_n = 256;
z_slice_n = 256;

phantom_full = zeros(x_slice_n,y_slice_n,z_slice_n);

switch category_phantom
    case 1
        size_phantom = [13 20 16];
        start_x = ceil(((x_slice_n*xyz_resolution/2) - (size_phantom(1)/2))./xyz_resolution);
        end_x = floor(((x_slice_n*xyz_resolution/2) + (size_phantom(1)/2))./xyz_resolution);
        range_x = start_x:end_x;

        start_y = ceil(((y_slice_n*xyz_resolution/2) - (size_phantom(2)/2))./xyz_resolution);
        end_y = floor(((y_slice_n*xyz_resolution/2) + (size_phantom(2)/2))./xyz_resolution);
        range_y = start_y:end_y;

        start_z = ceil(((z_slice_n*xyz_resolution/2) - (size_phantom(3)/2))./xyz_resolution);
        end_z = floor(((z_slice_n*xyz_resolution/2) + (size_phantom(3)/2))./xyz_resolution);
        range_z = start_z:end_z;

        phantom_full(range_x , range_y ,range_z ) = 7;
        figure,imshow3Dfull(phantom_full)

        save('phantom_flat_homo_A2_hybrid.mat','phantom_full')
end
