scale = [61 31];
n_laser = 49;
B = zeros(scale(1),scale(2),n_laser);
for i = 1:n_laser
    B(:,:,i) = imresize((reg_ex_im(:,:,i+1) -reg_ex_im(:,:,1)),scale);
end
min_Ex = min(B(:));
max_Ex = max(B(:));
figure,
for i=1:n_laser
    subplot(7,7,i),
    B_det = B(:,:,i)';
    imagesc(B_det);
    caxis([min_Ex max_Ex]) 
    title(['Ex map with laser #'  num2str(i)]);
end
figure,
for i=1:n_laser
    B_det = B(:,:,i)';
    B_line = reshape(B_det,[],1);
    subplot(1,2,1),plot(exp_FMT.det_Ex(:,i))
    subplot(1,2,2),plot(B_line)
    pause(0.5);
end