function BoxPhantom = drawFroPhantom_new(p_dim, p_dl, i_pos, i_dim, i_val,level)
    %% illustration
    %**********************************************************************
    % p_dim: box dimension [x, y, z]' /mm
    % p_dl: [dx,dy,dz]' /mm
    % i_pos: [px1 px2 px3 ...
    %         py1 py2 py3 ...
    %         pz1 pz2 pz3 ...]
    % i_dim: [x1 x2 x3 ...
    %         y1 y2 y3 ...
    %         z1 z2 z3 ...]
    % i_val: [val1 val2 val3... val_background]'
    %    
    %    ____________________________________
    %   |                                    |
    %   |    |--------|                      |
    %   |    |        |       val_background |
    %   |    |  val1  |                      |
    %   |    |        |                      |
    %   |    |--------|     |----|           | 
    %   |                   |val2|           |
    %   |                   |----|           |
    %   |                                    |
    %   |____________________________________|
    %**********************************************************************
    %%
    grd = floor(p_dim./p_dl)+1;
    % Gussassion back_ground
%     mu=[10 10 8];
%     sigma = [0.8,0.4,0.6; 0.4,0.9,0.5; 0.6,0.5,0.8];
%     [x y z]=meshgrid(linspace(1,31,31)',linspace(1,11,16)',linspace(1,11,11)');
%     X =[x(:) y(:) z(:)];
%  
%     back_ground=3000.*mvnpdf(X,mu,sigma);
%     index = 31*16*4+28;
%     %[m,index]=max(back_ground);
%     back_ground(index) = 200; back_ground(index+100) = 200;
%     % back_ground(index+31*16:index+31*16+2) = 100:20:140;
%     back_ground =  back_ground-min(back_ground).*ones(size(back_ground));
%     BoxPhantom = reshape(back_ground,31,16,11);
    
    
     back_ground = i_val(length(i_val)) * ones(grd');
%% 高斯密度背景分布
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% mu = 10.*ones(1,3);                            %该高斯分布的均值向量
% sigma = [0.1 0 0; 0 0.1 0; 0 0 0.1];                 %该高斯分布的协方差矩阵
% 
% % 使用 mvnrnd 函数得到符合该高斯布的样本点
% data = mvnrnd(mu,sigma,5000);
% x = data(:,1);                              %样本点横坐标
% y = data(:,2);                              %样本点纵坐标
% z = data(:,3); 
% % 第2步.
% % 使用 mvnpdf 函数得到每个点的概率密度
% f = mvnpdf(data,mu,sigma);
% [xq,yq,zq] = meshgrid(linspace(min(x),max(x),grd(2)),linspace(min(y),max(y),grd(1)),linspace(min(z),max(z),grd(3)));
% % 第3步.
% % 使用  griddata 和 surf 函数 对已经得到的样本点网格化，其中插值法我们选择使用 'v4'
% [value]=griddata(x,y,z,f,xq,yq,zq,'nearest'); %插值
% value(ceil(end/4):ceil(end/4)+ceil(end/2)-1,ceil(end/4):ceil(end/4)+ceil(end/2)-1,ceil(end/4):ceil(end/4)+ceil(end/2)-1) = 0.8*...
%     value(ceil(end/4):ceil(end/4)+ceil(end/2)-1,ceil(end/4):ceil(end/4)+ceil(end/2)-1,ceil(end/4):ceil(end/4)+ceil(end/2)-1);
% back_ground(1:ceil(end/2),:,:) = value(ceil(end/4):ceil(end/4)+ceil(end/2)-1,:,:);
% back_ground(ceil(end/2):end,:,:) = value(ceil(end/4):ceil(end/4)+end-ceil(end/2),:,:);
% [value_max, value_indmax] = max(value(:));
% % % 最大值坐标
% [i, j, w] = ind2sub(size(value), value_indmax);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%最原始数据
% back_ground(ceil(grd(1)/8)-2:ceil(grd(1)/2)+ceil(grd(1)/8)-3,ceil(grd(1)/8):ceil(grd(2)/2)+ceil(grd(1)/8)-1,1:grd(3)) = 0.6*value +...
%     back_ground(ceil(grd(1)/8)-2:ceil(grd(1)/2)+ceil(grd(1)/8)-3,ceil(grd(1)/8):ceil(grd(2)/2)+ceil(grd(1)/8)-1,1:grd(3))
% back_ground(ceil(grd(1)/8)-2:ceil(grd(1)/2)+ceil(grd(1)/8)-3,end - ceil(grd(2)/2)-1:end-2,1:grd(3)) = 0.6*value +...
%    back_ground(ceil(grd(1)/8)-2:ceil(grd(1)/2)+ceil(grd(1)/8)-3,end - ceil(grd(2)/2)-1:end-2,1:grd(3));
% back_ground(end - ceil(grd(1)/2):end-1,ceil(grd(1)/8):ceil(grd(2)/2)+ceil(grd(1)/8)-1,1:grd(3)) = 0.8*value +...
%     back_ground(end - ceil(grd(1)/2):end-1,ceil(grd(1)/8):ceil(grd(2)/2)+ceil(grd(1)/8)-1,1:grd(3));
% back_ground(end - ceil(grd(1)/2):end-1,end - ceil(grd(2)/2)-1:end-2,1:grd(3)) = value+...
%     back_ground(end - ceil(grd(1)/2):end-1,end - ceil(grd(2)/2)-1:end-2,1:grd(3));
% back_ground(ceil(grd(1)/4):ceil(grd(1)/2)+ceil(grd(1)/4)-1,ceil(grd(2)/4):ceil(grd(2)/2)+ceil(grd(2)/4)-1,1:grd(3)) = value+...
%     back_ground(ceil(grd(1)/4):ceil(grd(1)/2)+ceil(grd(1)/4)-1,ceil(grd(2)/4):ceil(grd(2)/2)+ceil(grd(2)/4)-1,1:grd(3));
% back_ground = back_ground + min(back_ground(:))*ones(grd');
%%%%%%%%%%%%%%%%%%%%
% value2 = value1(1:grd(1),1:grd(2),1:grd(3));
% value3 = value1(end-grd(1):end-1,end-grd(2):end-1,end-grd(3):end-1);
% [value2_max, value2_indmax] = max(value2(:));

% % 最大值坐标
% [i, j, w] = ind2sub(size(value2), value2_indmax);
% value2(i-1:i+1, j-1:j+1, w-1:w+1) = 1.2*value2(i-1:i+1, j-1:j+1, w-1:w+1);
% value2(i, j, w) = 1.2*value2(i, j, w);
% [value3_max, value3_indmax] = max(value3(:));
% % 最大值坐标
% [i, j, w] = ind2sub(size(value3), value3_indmax);
% value3(i-1:i+1, j-1:j+1, w-1:w+1) = 1.2*value3(i-1:i+1, j-1:j+1, w-1:w+1);
% value3(i, j, w) = 1.2*value3(i, j, w);
% value = value2 + value3 ;
% % 最大值统计
% [value1_max, value1_indmax] = max(value1(:))
% % 最大值坐标
% [i, j, w] = ind2sub(size(value1), value1_indmax);
% value1(i, j, w) = 2*value1(i, j, w);


%% 几个稀疏的点

% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % for i = 1:grd(3)
% %     temp = min(i,grd(3)-i)*10.*Z;
% %     back_ground(:,:,i) = temp-min(temp(:)).*ones(grd(1),grd(2));
% % end
% % temp = back_ground(:,:,6);
% % temp = temp(:);
% % [C,index] = max(temp);
% % temp(index) = C*20;
% % back_ground(:,:,6) = reshape(temp,grd(1),grd(2));
% % back_ground = zeros(grd');
% back_ground(8,12,3) = 100;
% 
% % back_ground(24,8,3) = 100; 
% 
% back_ground(16,4,3) = 100;
% %back_ground(15,10,2) = 100;
% BoxPhantom = back_ground;


%% fluo_bar
%      BoxPhantom = back_ground;
%     if isempty(i_pos)==0 && isempty(i_dim)==0
%         in_grd_1 = floor((i_pos-i_dim/2)./repmat(p_dl,1,size(i_pos,2)))+1;
%         in_grd_2 = floor((i_pos+i_dim/2)./repmat(p_dl,1,size(i_pos,2)))+1;
%         % error correction
%         for i =1:3
%            in_grd_1(i,in_grd_1(i,:)<1)=1;
%            in_grd_2(i,in_grd_2(i,:)<1)=1;
%            in_grd_1(i,in_grd_1(i,:)>grd(i))=grd(i);
%            in_grd_2(i,in_grd_2(i,:)>grd(i))=grd(i);
%         end
%         for i =1:size(i_pos,2)
%             BoxPhantom(in_grd_1(1,i):in_grd_2(1,i)-1,in_grd_1(2,i):in_grd_2(2,i),...
%                 in_grd_1(3,i):in_grd_2(3,i)) = i_val(i);
%         end
%     end

%% phantom 1
% back_ground(ceil(grd(1)/2),ceil(grd(2)/2)-ceil(grd(2)/4):ceil(grd(2)/2)+ceil(grd(2)/4),end-3) = 100;
% BoxPhantom = back_ground;
%% simulation 

    
% % diffuse scene
% load x_display.mat;
% % x_display(1:ceil(3*end/5),:,:) = 0.7*x_display(1:ceil(3*end/5),:,:);
% x_display(:,:,2) = x_display(:,:,3)/2;
% x_display(:,:,1) = x_display(:,:,2)/2;
% x_display(:,:,end) = x_display(:,:,end-1)/2;
% maxvalue = max(x_display(:));
% x_display(x_display>(maxvalue/2)) = 2/3*x_display(x_display>(maxvalue/2));
% back_ground = x_display;
% BoxPhantom = back_ground;
% BoxPhantom(28,9,7) = 10;
% BoxPhantom(15,14,5) = 10;

%% 同质的花生模型 
BoxPhantom = back_ground;
fltgt_g = zeros(grd');
nblobs = 2;
fcontrast = [0.1 0.1];
blbcx = grd(1)*[0.32 0.68]; % xcentre of blobs 0-41
blbcy = grd(2)*[0.5,0.5]; % ycentre of blobs 0-21
blbcz = grd(3)*[0.5,0.5]; % zcentre of blobs 0-11
blbrad = 2*grd(3)*[0.35 0.35]; % radius of blobs
% blbrad = 2*grd(3)*[0.3 0.3]; % radius of blobs
for i = 1:grd(1)
  for j = 1:grd(2)
      for m = 1:grd(3)
          for k = 1:nblobs
            if( (i-blbcx(k))^2 + (j-blbcy(k))^2 + 4*(m-blbcz(k))^2 < blbrad(k)^2)
                fltgt_g(i,j,m) = fcontrast(k);
            end
          end
      end
  end
end
back_ground = fltgt_g;
% 
BoxPhantom = back_ground;
% BoxPhantom(15,30,9) = 8;
% BoxPhantom(35,25,5) = 8;
% BoxPhantom(ceil((grd(2)-1)*2/3):ceil((grd(2)+2)*2/3)-1,ceil((grd(1)-1)/2):ceil((grd(1)+1)/2),ceil((grd(3)-1)/2):ceil((grd(3)+1)/2)-1) = 2;
BoxPhantom(37:38 ,27, 8) = 8;

end