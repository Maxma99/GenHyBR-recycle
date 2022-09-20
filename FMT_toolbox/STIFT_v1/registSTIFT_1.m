% v1
function obj = registSTIFT_1(obj)
    meshing_im_str = load([obj.project_directory '\mesh_top.mat']);
    fn = fieldnames(meshing_im_str);
    mesh_im = meshing_im_str.(fn{1});    
%     wt_im = imread([obj.measurement.wt_folder obj.measurement.wt_file]);

%% whitelight
testSequence = '.\testData18';
wt_im = imread([testSequence,'\whitelight_undist.png']);
load([testSequence,'\T_p2c.mat']);
    wt_im = double(wt_im);
%     wt_im_show = wt_im./max(wt_im(:));
%     wt_im_show = histeq(wt_im_show);
    %h=cpselect(wt_im_show,mesh_im);
    D = size(mesh_im);
    fixedPoints0 = [0.001,0.001,1;0.001,0.999,1;0.999,0.001,1;0.999,0.999,1;];
    movingPoints0 = fixedPoints0*T_p2c;
    movingPoints = movingPoints0(:,1:2);
    fixedPoints = size(mesh_im,1)*fixedPoints0(:,1:2);
%     [movingPoints,fixedPoints] = cpselect(wt_im_show,mesh_im,'Wait',true);
    t = cp2tform(movingPoints,fixedPoints,'similarity');
    
    % transform whitelight image & save
    reg_wt_im = imtransform(wt_im,t,'XData',[1 D(2)],'YData',[1 D(1)]);
    %figure,imagesc(reg_wt_im);
    %title('registered white light image')
    %axis image  
    save([obj.project_directory '\reg_wt_im.mat'],'reg_wt_im');
 
 %% excitation
    % transform the excitation map & save
    rawDataPath = [testSequence,'\RawData_640_undistortion\'];
x_range = [0.15,0.85];
y_range = [0.15,0.85];
nx = 10;
ny = 10;
[det_x,det_z] = meshgrid(linspace(x_range(1),x_range(2),nx),linspace(y_range(1),y_range(2),ny));
scan_points_xy = [reshape(det_x,[],1) reshape(det_z,[],1)];
scan_points_xy_1 = [scan_points_xy,0.01*ones(size(scan_points_xy,1),1)];

% load([testSequence,'\T_p2c_undistortion.mat']);
% load([testSequence,'\T_p2c.mat']);
% scan_points_camera_undistortion = scan_points_xy_1* T_p2c;
% load([testSequence,'\cameraParams.mat']);
count = size(scan_points_xy,1);
reg_ex_im = zeros([D count]);
% ambient_img = imread([testSequence,'\ambientlight_undistortion.png']);
ambient_img = imread([testSequence,'\ambientlight_640_undistortion.png']);
tic
for i=1:count
% imgName = [rawDataPath,num2str(scan_points_xy(i,1)), ',', num2str(scan_points_xy(i,2)),'_undistortion.png'];
imgName = [rawDataPath,num2str(scan_points_xy(i,2)), ',', num2str(scan_points_xy(i,1)),'_undistortion.png'];
tmp = imread(imgName);
% [tmp_undst,newOrigin] = undistortImage(tmp,cameraParams);
tmp_undst = tmp;
temp_img = tmp_undst - ambient_img;
reg_ex_im(:,:,i)= imtransform(temp_img,t,'XData',[1 D(2)],'YData',[1 D(1)]);
end
toc
    %figure,imshow3Dfull(trans_ex_im);
    save([obj.project_directory '\reg_ex_im.mat'],'reg_ex_im');
 %% emission   
    % transform the emission map & save
    rawDataPath = [testSequence,'\RawData_680_undistortion\'];
    ambient_img = imread([testSequence,'\ambientlight_680_undistortion.png']);
    reg_em_im = zeros([D count]);
    tic
for i=1:count
% imgName = [rawDataPath,num2str(scan_points_xy(i,1)), ',', num2str(scan_points_xy(i,2)),'_undistortion.png'];
imgName = [rawDataPath,num2str(scan_points_xy(i,2)), ',', num2str(scan_points_xy(i,1)),'_undistortion.png'];
tmp = imread(imgName);
% [tmp_undst,newOrigin] = undistortImage(tmp,cameraParams);
tmp_undst = tmp;
temp_img = tmp_undst - ambient_img;
reg_em_im(:,:,i)= imtransform(temp_img,t,'XData',[1 D(2)],'YData',[1 D(1)]);
end
toc
    %figure,imshow3Dfull(trans_ex_im);
    save([obj.project_directory '\reg_em_im.mat'],'reg_em_im');

end
