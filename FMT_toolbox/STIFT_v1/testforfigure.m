%%
% load('D:\1 code\01_project_s_decom\STIFT_v1\0_record\phantom_3(testData16)\data\l1FISTA.mat')
% load('D:\1 code\01_project_s_decom\STIFT_v1\0_record\phantom_3(testData16)\data\l1l2.mat')
% load('D:\1 code\01_project_s_decom\STIFT_v1\0_record\phantom_3(testData16)\data\ground_truth.mat')
% load('D:\1 code\01_project_s_decom\STIFT_v1\0_record\phantom_3(testData16)\data\l2CG.mat')
% load('D:\1 code\01_project_s_decom\STIFT_v1\0_record\phantom_3(testData16)\data\s1.mat')
% load('D:\1 code\01_project_s_decom\STIFT_v1\0_record\phantom_3(testData16)\data\s2.mat')
% load('D:\1 code\01_project_s_decom\STIFT_v1\0_record\phantom_3(testData16)\data\s.mat')
%%
s = x_out;
gt = fluo_g_disp;

% l1FISTA = l1FISTA2;
l1FISTA(l1FISTA<0) = 0;
l1l2(l1l2<0) = 0;
l2CG(l2CG<0) = 0;
s(s<0) = 0;
s1(s1<0) = 0;
s2(s2<0) = 0;

% l1FISTA = expand(l1FISTA);
% l1l2 = expand(l1l2);
% l2CG = expand(l2CG);
% s = expand(s);
% s1 = expand(s1);
% s2 = expand(s2);
% gt = fluo_g_disp;
index = [37,27,8];
gt_line = gt(:,index(2),index(3));
l1F_line = l1FISTA(:,index(2),index(3));
l1l2_line = l1l2(:,index(2),index(3));
l2CG_line = l2CG(:,index(2),index(3));
s_line = s(:,index(2),index(3));

%%
gt_line = gt_line/max(gt_line);
l1F_line = l1F_line/max(l1F_line);
l1l2_line = l1l2_line/max(l1l2_line);
l2CG_line = l2CG_line/max(l2CG_line);
s_line = s_line/max(s_line);
%%
figure
plot(gt_line,'linewidth',1.2,'color','black');
hold on
plot(l2CG_line,'-.','linewidth',1.2,'color','#9B30FF');
plot(l1F_line,'--','linewidth',1.2,'color','blue');
plot(l1l2_line,'--*','linewidth',1.2,'color','#00CDCD');
plot(s_line,'-p','linewidth',1.2,'color','#cd3700');
hold off
legend('Ground truth','LSCG($\ell_2$)','FISTA($\ell_1$)','MM($\ell_1 \& \ell_2$ )','SDHPM','interpreter','latex');
xlabel('Pixel')
ylabel('Norm.int(a.u.)')

%% cnr
%normalize
gt = gt/max(gt(:));
l1FISTA = l1FISTA/max(l1FISTA(:));
l1l2 = l1l2/max(l1l2(:));
l2CG = l2CG/max(l2CG(:));
s = s/max(s(:));
s1 = s1/max(s1(:));
s2 = s2/max(s2(:));

thr = 0.5;
l1F_cnr = computeCNR(l1FISTA,gt,thr);
l1l2_cnr = computeCNR(l1l2,gt,thr);
l2CG_cnr = computeCNR(l2CG,gt,thr);
s_cnr = computeCNR(s,gt,thr);

% cnr = [l2CG_cnr,l1F_cnr,l1l2_cnr,s_cnr];
cnr = [l1l2_cnr,l1F_cnr,l2CG_cnr,s_cnr];
cnr = cnr/max(cnr);

%% mse

l1F_rmse = sqrt(mse(gt(:),l1FISTA(:)));
l1l2_rmse = sqrt(mse(gt(:),l1l2(:)));
l2CG_rmse = sqrt(mse(gt(:),l2CG(:)));
s_rmse = sqrt(mse(gt(:),s(:)));
% rmse=[l2CG_rmse,l1F_rmse,l1l2_rmse,s_rmse];
rmse=[l1l2_rmse,l1F_rmse,l2CG_rmse,s_rmse];
rmse = rmse/max(rmse);

%% dice
% gt(37:39,27:29,7:9)=4;
s1_gt = gt; 
s1_gt(37:38 ,27, 8) = s1_gt(36 ,27, 8); % delete sparse objects
s1_l1F_dice = computeDice(l1FISTA,s1_gt,thr);
s1_l1l2_dice = computeDice(l1l2,s1_gt,thr);
s1_l2CG_dice = computeDice(l2CG,s1_gt,thr);
s1_dice = computeDice(s1,s1_gt,thr);
% dice1 = [s1_l2CG_dice,s1_l1F_dice,s1_l1l2_dice,s1_dice];
dice1 = [s1_l1l2_dice,s1_l1F_dice,s1_l2CG_dice,s1_dice];
% dice = [0.0158    0.0582    0.0418    0.1569];
dice1 = dice1/max(dice1);

s2_gt = zeros(size(gt));
s2_gt(37:38 ,27, 8) = 1;
s2_l1F_dice = computeDice(l1FISTA,s2_gt,thr);
s2_l1l2_dice = computeDice(l1l2,s2_gt,thr);
s2_l2CG_dice = computeDice(l2CG,s2_gt,thr);
s2_dice = computeDice(s2,s2_gt,thr);
% dice2 = [s2_l2CG_dice,s2_l1F_dice,s2_l1l2_dice,s2_dice];
dice2 = [s2_l1l2_dice,s2_l1F_dice,s2_l2CG_dice,s2_dice];
% dice = [0.0158    0.0582    0.0418    0.1569];
dice2 = dice2/max(dice2);

figure
y=[rmse;cnr;dice1;dice2];
b=bar(y);
grid on;
ch = get(b,'children');
set(gca,'XTickLabel',{'RMSE','CNR','Dice(smooth region)','Dice(sparse object)'})
b(1).FaceColor  = '#9B30FF';
b(2).FaceColor  = 'blue';
b(3).FaceColor  = '#00CDCD';
b(4).FaceColor  = '#cd3700';
legend('LSCG($\ell_2$)','FISTA($\ell_1$)','MM($\ell_1 \& \ell_2$ )','SDHPM','interpreter','latex');
ylabel('Performance (normalized)');


