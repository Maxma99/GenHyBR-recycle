% v1
function obj = loadMeasurement_1(obj)
    tic
    load_Ex = 1;
    fprintf(1,'loading Ex/Em from raw data...\n');
    scale = obj.detector.d_number;
    n_laser = size(obj.laser.l_position,1);
    if load_Ex == 1
        %% read Ex images
        B = zeros(scale(1),scale(2),n_laser);
        temp_load = load([obj.project_directory '\reg_ex_im.mat']);
        fn = fieldnames(temp_load);
        reg_ex_im = temp_load.(fn{1});
        for i = 1:n_laser
            B(:,:,i) = imresize((reg_ex_im(:,:,i+1) -reg_ex_im(:,:,1)),scale);
            B_det = B(:,:,i)';
            B_line = reshape(B_det,[],1);
            obj.det_Ex(:,i) = B_line;
        end
    end
    %% read Em images
    B = zeros(scale(1),scale(2),n_laser);
    temp_load = load([obj.project_directory '\reg_em_im.mat']);
    fn = fieldnames(temp_load);
    reg_em_im = temp_load.(fn{1});
    for i = 1:n_laser
        B(:,:,i) = imresize((reg_em_im(:,:,i+1) -reg_em_im(:,:,1)),scale);
        B_det = B(:,:,i)';
        B_line = reshape(B_det,[],1);
        obj.det_Em(:,i) = B_line;
    end
    fprintf(1,'done\n');
    toc
end