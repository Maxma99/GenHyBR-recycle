function obj = drawFluoSettingSTIFT_1(obj)
    if isempty(obj.figHdl_FluoSetting)
        obj.figHdl_FluoSetting = figure;
    else
        if ~isfield(obj.figHdl_FluoSetting,'number')
            clear obj.figHdl_FluoSetting
            obj.figHdl_FluoSetting = figure;
        else
            figure(obj.figHdl_FluoSetting);
            clf(obj.figHdl_FluoSetting,'reset');
        end
    end
    draw_fluo_method = 2;
    switch draw_fluo_method
        case 1 % translate from mesh
            Grid_Mesh_basis = toastBasis(obj.toast_mesh, obj.grd);
            fluo_g = Grid_Mesh_basis.Map('M->S', obj.fluo.fluoDis_m);
            fluo_g_disp=reshape(fluo_g,obj.grd);
            imshow3Dfull(fluo_g_disp);
            colormap(hot);
            title('ground truth fluo')
            xlabel('x/mm');
            ylabel('y/mm');
            zlabel('z/mm');
        case 2 % directly from grid
            fluo_g_disp = drawFroPhantom_new(obj.phantom_parameter.dim, ...
            obj.phantom_parameter.dl, ...
            obj.fluo.pos, ...
            obj.fluo.dim, ...
            obj.fluo.list_con,1000);
        save([obj.data_buffer_directory '/ground_truth.mat'],'fluo_g_disp','-v7.3');
        %obj.weighting_Matrix = [obj.data_buffer_directory '/ground_truth.mat']; 
            imshow3Dfull(fluo_g_disp);
            %colormap(hot);
            title('ground truth fluo')
            xlabel('x/mm');
            ylabel('y/mm');
            zlabel('z/mm');          
    end
end

% for i=1:330
%     for j = 1:330
%         for k = 1:90
%             Vq(i,j,k) = fluo_g_disp(floor((i+5)/6),floor((j+5)/6),floor((k+5)/6));
%         end
%     end
% end
