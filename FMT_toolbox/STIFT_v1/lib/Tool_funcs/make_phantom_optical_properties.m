%% Pre-processing for three phantoms
% Wuwei Ren, 2017 Jun 22
% category 1: homogeneous slab phanom 
% category 2: inhomogeneous slab phanom 
% category 3: inhomogeneous semicylindrical phanom 
close all;
clear all;
category_phantom = 3;
xyz_resolution = 20/128; % MRI xy resolution,mm
x_slice_n = 256;
y_slice_n = 256;
z_slice_n = 400;

phantom_full = zeros(x_slice_n,y_slice_n,z_slice_n);

switch category_phantom
    case 1
        size_phantom = [15 30 60];
        start_x = ceil(((x_slice_n*xyz_resolution/2) - (size_phantom(1)/2))./xyz_resolution);
        end_x = floor(((x_slice_n*xyz_resolution/2) + (size_phantom(1)/2))./xyz_resolution);
        range_x = start_x:end_x;

        start_y = ceil(((y_slice_n*xyz_resolution/2) - (size_phantom(2)/2))./xyz_resolution);
        end_y = floor(((y_slice_n*xyz_resolution/2) + (size_phantom(2)/2))./xyz_resolution);
        range_y = start_y:end_y;

        start_z = ceil(((z_slice_n*xyz_resolution/2) - (size_phantom(3)/2))./xyz_resolution);
        end_z = floor(((z_slice_n*xyz_resolution/2) + (size_phantom(3)/2))./xyz_resolution);
        range_z = start_z:end_z;

        phantom_full(range_x , range_y ,range_z ) = 5;
        figure,imshow3Dfull(phantom_full)

        save('phantom_flat_homo.mat','phantom_full')
    case 2
        size_phantom = [15 30 60];
        start_x = ceil(((x_slice_n*xyz_resolution/2) - (size_phantom(1)/2))./xyz_resolution);
        mid_x = round((x_slice_n*xyz_resolution/2)./xyz_resolution);
        end_x = floor(((x_slice_n*xyz_resolution/2) + (size_phantom(1)/2))./xyz_resolution);
        range_x_1 = start_x:mid_x;
        range_x_2 = mid_x:end_x;

        start_y = ceil(((y_slice_n*xyz_resolution/2) - (size_phantom(2)/2))./xyz_resolution);
        end_y = floor(((y_slice_n*xyz_resolution/2) + (size_phantom(2)/2))./xyz_resolution);
        range_y = start_y:end_y;

        start_z = ceil(((z_slice_n*xyz_resolution/2) - (size_phantom(3)/2))./xyz_resolution);
        end_z = floor(((z_slice_n*xyz_resolution/2) + (size_phantom(3)/2))./xyz_resolution);
        range_z = start_z:end_z;

        phantom_full(range_x_1 , range_y ,range_z ) = 6;
        phantom_full(range_x_2 , range_y ,range_z ) = 5;       
        figure,imshow3Dfull(phantom_full)

        save('phantom_flat_inhomo.mat','phantom_full')
    case 3
        size_phantom = [15 30 60];
        radius_phantom_2 = round(size_phantom(1)*size_phantom(1)/(xyz_resolution*xyz_resolution));
        
        start_x = ceil(((x_slice_n*xyz_resolution/2) - (size_phantom(1)/2))./xyz_resolution);
        mid_x = round((x_slice_n*xyz_resolution/2)./xyz_resolution);
        end_x = floor(((x_slice_n*xyz_resolution/2) + (size_phantom(1)/2))./xyz_resolution);
        range_x = start_x:end_x;
 
        start_y = ceil(((y_slice_n*xyz_resolution/2) - (size_phantom(2)/2))./xyz_resolution);
        end_y = floor(((y_slice_n*xyz_resolution/2) + (size_phantom(2)/2))./xyz_resolution);
        range_y = start_y:end_y;
        mid_y = round((y_slice_n*xyz_resolution/2)./xyz_resolution);
        
        start_z = ceil(((z_slice_n*xyz_resolution/2) - (size_phantom(3)/2))./xyz_resolution);
        end_z = floor(((z_slice_n*xyz_resolution/2) + (size_phantom(3)/2))./xyz_resolution);
        range_z = start_z:end_z;
        
        temp_img = zeros(x_slice_n,y_slice_n);
        for k = range_z
            temp_img = zeros(x_slice_n,y_slice_n);
            for j = range_y
                for i = range_x
                    dis_to_center_2 = (i - end_x)^2 + (j - mid_y)^2;
                    if (dis_to_center_2<radius_phantom_2)&(i<=end_x) 
                        if i>=mid_x
                            temp_img(i,j) = 5;
                        else
                            temp_img(i,j) = 6;
                        end
                    end
                end
            end
            phantom_full(:,:,k) = temp_img;
        end     
        figure,imshow3Dfull(phantom_full)
        save('phantom_curved_inhomo.mat','phantom_full')
end
