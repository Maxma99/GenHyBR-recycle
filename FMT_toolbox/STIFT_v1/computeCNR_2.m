function cnr = computeCNR_2(volume,gt,thr)
grd = size(gt);
cnt_roi =0;
cnt_rob =0;
if thr == 1
    for i=1:grd(1)
        for j=1:grd(2)
            for k = 1:grd(3)
                if i>18 && i<22 && j>26 && j<30 && k>6 && k<10
                    cnt_roi = cnt_roi +1;
                    ROI(cnt_roi) = volume(i,j,k);
                else
                    cnt_rob = cnt_rob + 1;
                    ROB(cnt_rob) = volume(i,j,k);
                end
            end
        end
    end
else
    for i=1:grd(1)
        for j=1:grd(2)
            for k = 1:grd(3)
                if gt(i,j,k)==0
                    cnt_rob = cnt_rob + 1;
                    ROB(cnt_rob) = volume(i,j,k);
                else
                    cnt_roi = cnt_roi + 1;
                    ROI(cnt_roi) = volume(i,j,k);
                end
            end
        end
    end
end

mu_ROI = mean(ROI);
mu_ROB = mean(ROB);
var_ROI = var(ROI);
var_ROB = var(ROB);

cnr = abs(mu_ROI - mu_ROB)/sqrt(var_ROI+var_ROB);
%     temp = zeros(n,1);
%     for i=1:total
%         if(image(i)>thr)
%             ROI(j) = image(i);
%             temp(j) = gt(i);
%             j = j+1;
%         end
%     end
%     cnr = 10*log10(norm(temp)^2/norm(temp-ROI)^2);
end