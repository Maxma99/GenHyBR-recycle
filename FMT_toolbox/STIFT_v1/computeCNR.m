function cnr = computeCNR(volume,gt,thr)
    image = volume(:);
    gt = gt(:);
    maxvalue = max(image);
    total = length(image);
    n = sum(image>thr);
    ROI = zeros(n,1);
    ROB = zeros(total-n,1);
    j = 1;
    k = 1;
    for i=1:total
        if(image(i)>thr)
            ROI(j) = image(i);
            j = j+1;
        else
            ROB(k) = image(i);
            k = k+1;
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