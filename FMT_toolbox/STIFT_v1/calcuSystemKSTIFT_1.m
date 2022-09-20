% v1
function obj = calcuSystemKSTIFT_1(obj)
    obj.calcuMvecSTIFT;
    obj.calcuQvecSTIFT;
    tic
    fprintf(1,'Generating a system matrix K from TOAST++...\n');
    system_K = dotSysmat (obj.toast_mesh,obj.optiProp.o_mua,obj.optiProp.o_mus,...
       obj.optiProp.o_ref,obj.laser.l_freq); % system matrix generation
    system_K = real(system_K);
    save([obj.data_buffer_directory '\system_K.mat'],'system_K');
    obj.system_K = [obj.data_buffer_directory '\system_K.mat'];
    fprintf(1,'done\n');
    toc
end