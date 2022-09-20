function dice = computeDice(volume,gt,thr)
% thr volume： normalized data
thr_gt = thr * max(gt(:));
thr_volume = thr * max(volume(:));
grd = size(gt);
mask_gt = gt(:)>thr_gt;
mask_volume = volume(:)>thr_volume;
cnt = 0;
for i = 1:prod(grd)
    if mask_gt(i)==1 && mask_volume(i)==1
        cnt = cnt + 1;
    end
end
dice = 2*cnt/(sum(mask_gt)+sum(mask_volume));
end