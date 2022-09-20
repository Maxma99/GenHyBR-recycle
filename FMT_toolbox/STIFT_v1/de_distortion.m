% 根据cameraParams 参数进行去畸变

%% calibraion(相机振镜)图像去畸变
testSequence = '.\testData14';
load([testSequence,'\cameraParams.mat'])
readPath = [testSequence,'\calibrationImgs\'];
savePath = [testSequence,'\calibrationImgs_undistortion\'];
if ~exist(savePath,'dir')
	mkdir(savePath)
end
delete([savePath, '*.png']);

files = dir([readPath,'*.png']);
for i = 1:length(files)
    imgName = files(i).name;
    newName = [extractBefore(imgName,'.png'),'_undistortion.png'];
    img = imread([readPath,imgName]);
    [J,newOrigin] = undistortImage(img,cameraParams);
    imwrite(J, [savePath,newName]);
end




%% raw data de-distortion
testSequence = '.\testData18';
load([testSequence,'\cameraParams.mat'])
readPath = [testSequence,'\RawData_640\'];
savePath = [testSequence,'\RawData_640_undistortion\'];
if ~exist(savePath,'dir')
	mkdir(savePath)
end
delete([savePath, '*.png']);

files = dir([readPath,'*.png']);
for i = 1:length(files)
    imgName = files(i).name;
    newName = [extractBefore(imgName,'.png'),'_undistortion.png'];
    img = imread([readPath,imgName]);
    [J,newOrigin] = undistortImage(img,cameraParams);
    imwrite(J, [savePath,newName]);
end

%% while light / ambient light image de-distortion
img = imread('ambientlight_640.png');
[J,newOrigin] = undistortImage(img,cameraParams);
imwrite(J,'ambientlight_640_undistortion.png');
img = imread('ambientlight_680.png');
[J,newOrigin] = undistortImage(img,cameraParams);
imwrite(J,'ambientlight_680_undistortion.png');
% figure()
% imshow(img)
% figure()
% imshow(J)
% kk = abs(double(img)-double(J));