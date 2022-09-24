% v1
function obj = geneWMSTIFT_1(obj,system_K_file)
%%
    tic
    fprintf(1,'Constructing the weighting matrix...\n');
%% load files
    %% load system_K
    system_K_s = load(system_K_file);
    system_K = system_K_s.system_K;
    clear system_K_s
    %%%%%
    %% load mvec
    mvec_s = load(obj.mvec);
    mvec = mvec_s.mvec;
    clear mvec_s
    %% load qvec
    qvec_s = load(obj.qvec);
    qvec = qvec_s.qvec;
    clear qvec_s
%% calculate the weighting matrix
    lgamma = reshape(full(obj.det_Ex),[],1);
    xdata = real(lgamma(obj.mea_mask));
    asol1 = system_K\mvec; % adjoint matrix
    nQ = size(qvec,2); % the number of sources
    nM = size(mvec,2); % the number of detectors
    %nsol = prod(obj.recon_grd); 
    Jnew = zeros(sum(obj.mea_mask),sum(obj.sol_mask));
    Grid_Mesh_basis = toastBasis(obj.toast_mesh, obj.recon_grd);
    k =1;
    for i = 1:nQ
      for j = 1:nM
         if obj.mea_mask(( i-1)*nM +j) 
            tmp_h = full(obj.phi_Ex(:,i).*asol1(:,j));    
            tmp_b = Grid_Mesh_basis.Map('M->B', tmp_h);
            %tmp_b = Grid_Mesh_basis.Map('M->S', tmp_h);
            Jnew(k,:) = tmp_b;
            k = k+1;
         end
      end
    end
    for i = 1: sum(obj.mea_mask)
    	Jnew(i,:) = Jnew(i,:)/xdata(i);
    end
    weighting_Matrix = Jnew;
    save([obj.data_buffer_directory '/weighting_Matrix.mat'],'weighting_Matrix','-v7.3');

    %save([obj.data_buffer_directory '\weighting_Matrix.mat'],'Jnew');
    obj.weighting_Matrix = [obj.data_buffer_directory '/weighting_Matrix.mat'];    

    fprintf(1,'done\n');
    toc    
end