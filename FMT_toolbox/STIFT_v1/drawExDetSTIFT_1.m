% v1
function obj=drawExDetSTIFT_1(obj)
    if isempty(obj.figHdl_Ex)
        obj.figHdl_Ex = figure;
    else
        if ~isfield(obj.figHdl_Ex,'number')
            clear obj.figHdl_Ex
            obj.figHdl_Ex = figure;
        else
            figure(obj.figHdl_Ex);
            clf(obj.figHdl_Ex,'reset');
        end
    end
    
    n_laser = length(obj.laser.l_position(:,1));
    min_Ex=full(min(min(obj.det_Ex)));
    max_Ex=full(max(max(obj.det_Ex)));
    if n_laser == prod(obj.laser.l_number)
        for i=1:n_laser
            subplot(obj.laser.l_number(1),obj.laser.l_number(2),i),
            temp_img = reshape(obj.det_Ex(:,i),obj.detector.d_number(2),obj.detector.d_number(1));
            imagesc(abs(temp_img));
            caxis([min_Ex max_Ex]) 
            title(['Ex map with laser #'  num2str(i)]);
        end
    else
        subplot_index_x = floor(sqrt(n_laser));
        subplot_index_y =ceil(n_laser/subplot_index_x);
        for i=1:n_laser
            subplot(subplot_index_x,subplot_index_y,i),
            temp_img = reshape(obj.det_Ex(:,i),obj.detector.d_number(2),obj.detector.d_number(1));
            imagesc(temp_img);
            caxis([min_Ex max_Ex]) 
            title(['Ex map with laser #'  num2str(i)]);
        end        
    end
end