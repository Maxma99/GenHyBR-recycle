% v1
function obj = geneMeaMaskSTIFT_1(obj)
    tic
    fprintf(1,'Constructing a measurement mask for reconstruction...\n');
    if obj.mea_mask_option
        lgamma = reshape(full(obj.det_Ex),[],1);
        obj.mea_mask = (lgamma>(obj.mea_mask_thr*max(lgamma))); % a simple threshold method
        fprintf(1,'done\n');          
    else
        lgamma = reshape(full(obj.det_Ex),[],1);
        obj.mea_mask = true(length(lgamma),1); % all ones        
        fprintf(1,'done(not applied)\n');           
    end
    toc   
end