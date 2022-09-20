% v1
function obj = geneDisposeMeaMaskSTIFT_1(obj)
    tic
    fprintf(1,'Constructing a measurement mask for reconstruction...\n');
    
    measure_array_s = load(obj.measure_array);
    measure_array = measure_array_s.measure_array;
    clear measure_array_s
    
    if obj.dispose_mea_mask_option
        lgamma = measure_array;
        obj.dispose_mea_mask = (lgamma>(obj.dispose_mea_mask_thr*max(lgamma))); % a simple threshold method
        % 随机取样部分被删除的测量值
        temp = rand(length(obj.dispose_mea_mask(obj.dispose_mea_mask==0)),1);
        temp(temp>0.8) = 1;
        temp(temp<=0.8) = 0;
        obj.dispose_mea_mask(obj.dispose_mea_mask==0) = temp;
        fprintf(1,'done\n');          
    else
        lgamma = measure_array;
        obj.dispose_mea_mask = true(length(lgamma),1); % all ones        
        fprintf(1,'done(not applied)\n');           
    end
    toc   
end