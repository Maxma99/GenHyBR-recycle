function obj = calcuQvecSTIFT_1(obj)  % v4.1
    %obj.qvec_method = 'gaussian-neumann';% options: 'collimated beam', 'gaussian-neumann','gaussian-point';
    tic
    fprintf(1,'Calculating Qvec...\n');
    switch obj.qvec_method
        case 'gaussian-neumann'
            qvec = obj.laser.l_power*obj.toast_mesh.Qvec('Neumann','Gaussian',obj.laser.l_width);
        case 'gaussian-point'
            qvec = obj.laser.l_power*obj.toast_mesh.Qvec('Neumann','Point',obj.laser.l_width);
        case 'collimated beam'
            % under construction
    end
    save([obj.data_buffer_directory '\qvec.mat'],'qvec');
    obj.qvec = [obj.data_buffer_directory '\qvec.mat'];
    fprintf(1,'done\n');
    toc
end