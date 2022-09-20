function obj = drawReconSTIFT_1(obj)
    iso_int = 0.4;
    if isempty(obj.figHdl_Recon)
        obj.figHdl_Recon = figure;
    else
        if ~isfield(obj.figHdl_Recon,'number')
            clear obj.figHdl_Recon
            obj.figHdl_Recon = figure;
        else
            figure(obj.figHdl_Recon);
            clf(obj.figHdl_Recon,'reset');
        end
    end
    
    x_display=reshape(obj.fluorecon,obj.recon_grd);
    save([obj.data_buffer_directory '/x_display.mat'],'x_display','-v7.3');
    imshow3Dfull(x_display);
    colormap(hot);
    title('reconstruction result')
    xlabel('x/mm');
    ylabel('y/mm');
    zlabel('z/mm');
    if ~isempty(obj.figHdl_ExpSetting)
        obj_bound = max(obj.node);
        x_display_T = permute(x_display,[2,1,3]);
        [m,n,p] = size(x_display_T);
        dl = [obj_bound(2)  obj_bound(1) obj_bound(3)]./[m-1 n-1 p-1];
        [X,Y,Z] = meshgrid(0:(n-1),0:(m-1),0:(p-1));
        X=X*dl(2);
        Y=Y*dl(1);
        Z=Z*dl(3);
        x_display_T(x_display_T<(max(x_display_T(:))*0.3))=0; %to reduce the artefact
        figure(obj.figHdl_ExpSetting),hold on
        patch(isosurface(X,Y,Z,x_display_T,max(x_display_T(:))*iso_int), ...
        'facecolor', '#00ffff', 'edgecolor', 'none');
        hold off
    end 
end