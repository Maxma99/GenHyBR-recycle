% v1
function obj = calcuMvecSTIFT_1(obj) % v4.1
    %obj.mvec_method = 'free-space';% options: 'free-space', 'gaussian', 'ortho';
    tic
    fprintf(1,'Calculating Mvec ');    
    switch obj.mvec_method
        case 'free-space'
            fprintf(1,'calculating travel matrix: Free-space propagation...\n');
            % to be constructed
        case 'gaussian' 
            fprintf(1,'with travel matrix: Gaussian...\n');
            mvec = obj.toast_mesh.Mvec('Gaussian', ...
                (obj.detector.d_dim(1)+obj.detector.d_dim(2))/2, obj.optiProp.o_ref);
        case 'ortho'
    end
    save([obj.data_buffer_directory '\mvec.mat'],'mvec');
    obj.mvec = [obj.data_buffer_directory '\mvec.mat'];
    fprintf(1,'done\n');
    toc
end