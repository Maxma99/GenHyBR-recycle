% v1
function obj = fluoSettingSTIFT_1(obj,Grid_Mesh_basis)
    %         fluo.fluoDis_m [real column vector]:
    %             nodal fluorescence dye distribution value
    %         fluo.method [string]:
    %             'regu' or 'load'
    switch obj.fluo.method
        case 'regu' % build the optical properties
            fluoDis_g = drawFroPhantom_new(obj.phantom_parameter.dim, ...
                obj.phantom_parameter.dl, ...
                obj.fluo.pos, ...
                obj.fluo.dim, ...
                obj.fluo.list_con,1000);
            obj.fluo.fluoDis_m  = Grid_Mesh_basis.Map ('B->M', ...
                reshape(fluoDis_g,[],1));% Basis Mapping
            %figure,imshow3Dfull(fluoDis_g);
        case 'load' % load the optical properties
    end
end   
