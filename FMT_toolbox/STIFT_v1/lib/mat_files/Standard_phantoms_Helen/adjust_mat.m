phantom_full = final_mri_256;
phantom_temp = permute(phantom_full,[3,2,1]);
phantom_temp_2 = flip(phantom_temp,3);
figure,imshow3Dfull(phantom_temp_2);
phantom_full= phantom_temp_2;