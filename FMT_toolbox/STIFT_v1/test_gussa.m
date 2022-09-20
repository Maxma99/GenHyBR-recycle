% clear;
% clc;
% 
% mu = zeros(1,2);                            %�ø�˹�ֲ��ľ�ֵ����
% sigma = [0.9 0.5; 0.5 0.9];                 %�ø�˹�ֲ���Э�������
% 
% % ��1��.
% % ʹ�� mvnrnd �����õ����ϸø�˹����������
% data = mvnrnd(mu,sigma,1000);
% x = data(:,1);                              %�����������
% y = data(:,2);                              %������������
% 
% % ��2��.
% % ʹ�� mvnpdf �����õ�ÿ����ĸ����ܶ�
% f = mvnpdf(data,mu,sigma);
% 
% % ��3��.
% % ʹ��  griddata �� surf ���� ���Ѿ��õ������������񻯣����в�ֵ������ѡ��ʹ�� 'v4'
% [X,Y,Z]=griddata(x,y,f,linspace(min(x),max(x),31)',linspace(min(y),max(y),16),'v4'); %��ֵ
% figure; surf(X,Y,Z); %������ 




 
%% By lyqmath @ Matlab������̳
% clc; clear all; close all;
% A = rand(5, 5, 5);
% % ���ֵͳ��
% [Amax, indmax] = max(A(:))
% % ���ֵ����
% [i, j, w] = ind2sub(size(A), indmax)
% % ��֤
% A(i, j, w)

%%
mu = 10.*ones(1,3);                            %�ø�˹�ֲ��ľ�ֵ����
sigma = [0.9 0 0; 0 0.9 0; 0 0 0.9];                 %�ø�˹�ֲ���Э�������

% ʹ�� mvnrnd �����õ����ϸø�˹����������
data = mvnrnd(mu,sigma,5000);
x = data(:,1);                              %�����������
y = data(:,2);                              %������������
z = data(:,3); 
% ��2��.
% ʹ�� mvnpdf �����õ�ÿ����ĸ����ܶ�
f = mvnpdf(data,mu,sigma);
[xq,yq,zq] = meshgrid(linspace(min(x),max(x),26),linspace(min(y),max(y),41),linspace(min(z),max(z),16));
% ��3��.
% ʹ��  griddata �� surf ���� ���Ѿ��õ������������񻯣����в�ֵ������ѡ��ʹ�� 'v4'
[value1]=griddata(x,y,z,f,xq,yq,zq,'nearest'); %��ֵ
