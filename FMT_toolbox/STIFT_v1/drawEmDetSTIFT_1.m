function obj=drawEmDetSTIFT_1(obj)
    if isempty(obj.figHdl_Em)
        obj.figHdl_Em = figure;
    else
        if ~isfield(obj.figHdl_Em,'number')
            clear obj.figHdl_Em
            obj.figHdl_Em = figure;
        else
            figure(obj.figHdl_Em);
            clf(obj.figHdl_Em,'reset');
        end
    end
    
    n_laser = length(obj.laser.l_position(:,1));
    min_Em=full(min(min(obj.det_Em)));
    max_Em=full(max(max(obj.det_Em)));
    if n_laser == prod(obj.laser.l_number)
        for i=1:n_laser
            subplot(obj.laser.l_number(1),obj.laser.l_number(2),i),
            temp_img = reshape(obj.det_Em(:,i),obj.detector.d_number(2),obj.detector.d_number(1));
            imagesc(temp_img);
            caxis([min_Em max_Em]); 
            title(['Em map with laser #'  num2str(i)]);
        end
        figure;
        %%%%%%%%%%%%%%%%%%% phantom 1 transmission
%     obj.det_Em(:,5) = [];
%     obj.det_Em(:,20) = [];
%     obj.det_Em(:,23) = [];
    %%%%%%%%%%%%%%%%%%%
%         for i=1:n_laser
%             subplot(obj.laser.l_number(1),obj.laser.l_number(2),i),
%             load('D:\1 code\STIFT_v1\pro\1_data\det_Em.mat')
%             temp_img = reshape(obj.det_Em(:,i),obj.detector.d_number(2),obj.detector.d_number(1));
%             temp_img = temp_img/max(temp_img(:));
%             temp_img_com = reshape(det_Em(:,i),obj.detector.d_number(2),obj.detector.d_number(1));
%             temp_img_com = temp_img_com/max(temp_img_com(:));
%             plot(temp_img(ceil(end/2),:));
%             hold on 
%             plot(temp_img_com(ceil(end/2),:));
%             hold off
%             legend('simulation','real data');
%         end
    else
        subplot_index_x = floor(sqrt(n_laser));
        subplot_index_y =ceil(n_laser/subplot_index_x);
        for i=1:n_laser
            subplot(subplot_index_x,subplot_index_y,i),
            temp_img = reshape(obj.det_Em(:,i),obj.detector.d_number(2),obj.detector.d_number(1));
            imagesc(temp_img);
            caxis([min_Em max_Em]); 
            title(['Em map with laser #'  num2str(i)]);
        end        
    end
end