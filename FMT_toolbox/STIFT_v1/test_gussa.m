% clear;
% clc;
% 
% mu = zeros(1,2);                            %该高斯分布的均值向量
% sigma = [0.9 0.5; 0.5 0.9];                 %该高斯分布的协方差矩阵
% 
% % 第1步.
% % 使用 mvnrnd 函数得到符合该高斯布的样本点
% data = mvnrnd(mu,sigma,1000);
% x = data(:,1);                              %样本点横坐标
% y = data(:,2);                              %样本点纵坐标
% 
% % 第2步.
% % 使用 mvnpdf 函数得到每个点的概率密度
% f = mvnpdf(data,mu,sigma);
% 
% % 第3步.
% % 使用  griddata 和 surf 函数 对已经得到的样本点网格化，其中插值法我们选择使用 'v4'
% [X,Y,Z]=griddata(x,y,f,linspace(min(x),max(x),31)',linspace(min(y),max(y),16),'v4'); %插值
% figure; surf(X,Y,Z); %画曲面 




 
%% By lyqmath @ Matlab中文论坛
% clc; clear all; close all;
% A = rand(5, 5, 5);
% % 最大值统计
% [Amax, indmax] = max(A(:))
% % 最大值坐标
% [i, j, w] = ind2sub(size(A), indmax)
% % 验证
% A(i, j, w)

%%
mu = 10.*ones(1,3);                            %该高斯分布的均值向量
sigma = [0.9 0 0; 0 0.9 0; 0 0 0.9];                 %该高斯分布的协方差矩阵

% 使用 mvnrnd 函数得到符合该高斯布的样本点
data = mvnrnd(mu,sigma,5000);
x = data(:,1);                              %样本点横坐标
y = data(:,2);                              %样本点纵坐标
z = data(:,3); 
% 第2步.
% 使用 mvnpdf 函数得到每个点的概率密度
f = mvnpdf(data,mu,sigma);
[xq,yq,zq] = meshgrid(linspace(min(x),max(x),26),linspace(min(y),max(y),41),linspace(min(z),max(z),16));
% 第3步.
% 使用  griddata 和 surf 函数 对已经得到的样本点网格化，其中插值法我们选择使用 'v4'
[value1]=griddata(x,y,z,f,xq,yq,zq,'nearest'); %插值
