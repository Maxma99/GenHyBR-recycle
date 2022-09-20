% v1
function obj = geneSolMaskSTIFT_1(obj)
    tic
    fprintf(1,'Constructing a solution mask for reconstruction...\n');
    if obj.sol_mask_option
        nsol = prod(obj.recon_grd); 
        obj.sol_mask = true(nsol,1); % to be constructed   
        fprintf(1,'done\n');          
    else
        nsol = prod(obj.recon_grd); 
        obj.sol_mask = true(nsol,1); % all ones        
        fprintf(1,'done(not applied)\n');           
    end
    toc   
end