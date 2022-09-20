% v1
function obj = drawExpSettingSTIFT_1(obj)
    if isempty(obj.figHdl_ExpSetting)
        obj.figHdl_ExpSetting = figure;
    else
        if ~isfield(obj.figHdl_ExpSetting,'number')
            clear obj.figHdl_ExpSetting
            obj.figHdl_ExpSetting = figure;
        else
            figure(obj.figHdl_ExpSetting);
            clf(obj.figHdl_ExpSetting,'reset');
        end
    end
    %% drawing mesh 
    plotmesh(obj.node(:,[1 2 3]),obj.face,'facealpha',0.3,'edgecolor','#636363');
%     title('FMT experimental setting')
    xlabel('x/mm')
    ylabel('y/mm')
    zlabel('z/mm')
    %% drawing laser
    hold on
    plot3(obj.laser.l_position(:,1),obj.laser.l_position(:,2),...
        obj.laser.l_position(:,3),'ro','MarkerFaceColor','r');
    % drawing detector
    if (obj.detector.d_normal_dA(1) <0)...
        &&(obj.detector.d_normal_dA(2)==0)...
        &&(obj.detector.d_normal_dA(3)==0)
        detector_orientation = 1 ;
    end
    if (obj.detector.d_normal_dA(2) <0)...
        &&(obj.detector.d_normal_dA(1)==0)...
        &&(obj.detector.d_normal_dA(3)==0)
        detector_orientation = 2 ;
    end
    if (obj.detector.d_normal_dA(3) <0)...
        &&(obj.detector.d_normal_dA(1)==0)...
        &&(obj.detector.d_normal_dA(2)==0)
        detector_orientation = 3 ;
    end
    if (obj.detector.d_normal_dA(1) >0)...
        &&(obj.detector.d_normal_dA(2)==0)...
        &&(obj.detector.d_normal_dA(3)==0)
        detector_orientation = 4 ;
    end
    if (obj.detector.d_normal_dA(2) >0)...
        &&(obj.detector.d_normal_dA(1)==0)...
        &&(obj.detector.d_normal_dA(3)==0)
        detector_orientation = 5 ;
    end
    if (obj.detector.d_normal_dA(3) >0)...
        &&(obj.detector.d_normal_dA(1)==0)...
        &&(obj.detector.d_normal_dA(2)==0)
        detector_orientation = 6 ;
    end
    switch detector_orientation
        case {1,4}
            for i = 1:length(obj.detector.d_virtualPos(:,1))
                d_x = [obj.detector.d_virtualPos(i,1) obj.detector.d_virtualPos(i,1)...
                    obj.detector.d_virtualPos(i,1) obj.detector.d_virtualPos(i,1)];
                
                d_y = [obj.detector.d_virtualPos(i,2)-obj.detector.d_dim(1)/2 ...
                    obj.detector.d_virtualPos(i,2)+obj.detector.d_dim(1)/2 ...
                    obj.detector.d_virtualPos(i,2)+obj.detector.d_dim(1)/2 ...
                    obj.detector.d_virtualPos(i,2)-obj.detector.d_dim(1)/2];
                
                d_z = [obj.detector.d_virtualPos(i,3)-obj.detector.d_dim(2)/2 ...
                    obj.detector.d_virtualPos(i,3)-obj.detector.d_dim(2)/2 ...
                    obj.detector.d_virtualPos(i,3)+obj.detector.d_dim(2)/2 ...
                    obj.detector.d_virtualPos(i,3)+obj.detector.d_dim(2)/2];
                col = 'b';
                patch(d_x,d_y,d_z,col);
            end           
        case {2,5}
             for i = 1:length(obj.detector.d_virtualPos(:,1))
                d_y = [obj.detector.d_virtualPos(i,2) obj.detector.d_virtualPos(i,2)...
                    obj.detector.d_virtualPos(i,2) obj.detector.d_virtualPos(i,2)];
                
                d_x = [obj.detector.d_virtualPos(i,1)-obj.detector.d_dim(1)/2 ...
                    obj.detector.d_virtualPos(i,1)+obj.detector.d_dim(1)/2 ...
                    obj.detector.d_virtualPos(i,1)+obj.detector.d_dim(1)/2 ...
                    obj.detector.d_virtualPos(i,1)-obj.detector.d_dim(1)/2];
                
                d_z = [obj.detector.d_virtualPos(i,3)-obj.detector.d_dim(2)/2 ...
                    obj.detector.d_virtualPos(i,3)-obj.detector.d_dim(2)/2 ...
                    obj.detector.d_virtualPos(i,3)+obj.detector.d_dim(2)/2 ...
                    obj.detector.d_virtualPos(i,3)+obj.detector.d_dim(2)/2];
                col = 'b';
                patch(d_x,d_y,d_z,col);
            end              
        case {3,6}
           for i = 1:length(obj.detector.d_virtualPos(:,1))
                d_z = [obj.detector.d_virtualPos(i,3) obj.detector.d_virtualPos(i,3)...
                    obj.detector.d_virtualPos(i,3) obj.detector.d_virtualPos(i,3)];
                
                d_x = [obj.detector.d_virtualPos(i,1)-obj.detector.d_dim(1)/2 ...
                    obj.detector.d_virtualPos(i,1)+obj.detector.d_dim(1)/2 ...
                    obj.detector.d_virtualPos(i,1)+obj.detector.d_dim(1)/2 ...
                    obj.detector.d_virtualPos(i,1)-obj.detector.d_dim(1)/2];
                
                d_y = [obj.detector.d_virtualPos(i,2)-obj.detector.d_dim(2)/2 ...
                    obj.detector.d_virtualPos(i,2)-obj.detector.d_dim(2)/2 ...
                    obj.detector.d_virtualPos(i,2)+obj.detector.d_dim(2)/2 ...
                    obj.detector.d_virtualPos(i,2)+obj.detector.d_dim(2)/2];
                col = 'b';
                patch(d_x,d_y,d_z,col);
            end    
    end
% 
    hold off
    % insert the whitelight image under construction
    %%%%%%%%%%%%%%
    hold on
    if obj.machine_mode == 3
        vari_wt = load([obj.project_directory '\reg_wt_im.mat']); % read the transformed ex image for recon recon recon!
        fn = fieldnames(vari_wt);
        wt_img = vari_wt.(fn{1});
        switch detector_orientation
            case 1        
                map_X = [min(obj.node(:,1)) min(obj.node(:,1)); min(obj.node(:,1)) min(obj.node(:,1))];
                map_Y = [min(obj.node(:,2)) max(obj.node(:,2)); min(obj.node(:,2)) max(obj.node(:,2))];
                map_Z = [min(obj.node(:,3)) min(obj.node(:,3)); max(obj.node(:,3)) max(obj.node(:,3))];              
            case 2
                map_X = [min(obj.node(:,1)) max(obj.node(:,1)); min(obj.node(:,1)) max(obj.node(:,1))];
                map_Y = [min(obj.node(:,2)) min(obj.node(:,2)); min(obj.node(:,2)) min(obj.node(:,2))];
                map_Z = [min(obj.node(:,3)) min(obj.node(:,3)); max(obj.node(:,3)) max(obj.node(:,3))];
            case 3
                map_X = [min(obj.node(:,1)) min(obj.node(:,1)); max(obj.node(:,1)) max(obj.node(:,1))];
                map_Y = [min(obj.node(:,2)) max(obj.node(:,2)); min(obj.node(:,2)) max(obj.node(:,2))];
                map_Z = [min(obj.node(:,3)) min(obj.node(:,3)); min(obj.node(:,3)) min(obj.node(:,3))];
            case 4
                map_X = [max(obj.node(:,1)) max(obj.node(:,1)); max(obj.node(:,1)) max(obj.node(:,1))];
                map_Y = [min(obj.node(:,2)) max(obj.node(:,2)); min(obj.node(:,2)) max(obj.node(:,2))];
                map_Z = [min(obj.node(:,3)) min(obj.node(:,3)); max(obj.node(:,3)) max(obj.node(:,3))];
            case 5
                map_X = [min(obj.node(:,1)) max(obj.node(:,1)); min(obj.node(:,1)) max(obj.node(:,1))];
                map_Y = [max(obj.node(:,2)) max(obj.node(:,2)); max(obj.node(:,2)) max(obj.node(:,2))];
                map_Z = [min(obj.node(:,3)) min(obj.node(:,3)); max(obj.node(:,3)) max(obj.node(:,3))];
            case 6
                map_X = [min(obj.node(:,1)) min(obj.node(:,1)); max(obj.node(:,1)) max(obj.node(:,1))];
                map_Y = [min(obj.node(:,2)) max(obj.node(:,2)); min(obj.node(:,2)) max(obj.node(:,2))];
                map_Z = [max(obj.node(:,3)) max(obj.node(:,3)); max(obj.node(:,3)) max(obj.node(:,3))];
        end
        surface('XData',map_X,...
        'YData',map_Y,...
        'ZData',map_Z,...
        'CData',wt_img',...
        'FaceColor','texturemap');
        colormap(gray)
        hold off
    end 
    
end