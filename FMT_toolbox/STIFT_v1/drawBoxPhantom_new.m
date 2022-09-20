function BoxPhantom = drawBoxPhantom_new(p_dim, p_dl, i_pos, i_dim, i_val)
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
    %BoxPhantom = max(randn(grd'),0);
    BoxPhantom = i_val(length(i_val)) * ones(grd');% background value
    if isempty(i_pos)==0 && isempty(i_dim)==0
        in_grd_1 = floor((i_pos-i_dim/2)./repmat(p_dl,1,size(i_pos,2)))+1;
        in_grd_2 = floor((i_pos+i_dim/2)./repmat(p_dl,1,size(i_pos,2)))+1;
        % error correction
        for i =1:3
           in_grd_1(i,in_grd_1(i,:)<1)=1;
           in_grd_2(i,in_grd_2(i,:)<1)=1;
           in_grd_1(i,in_grd_1(i,:)>grd(i))=grd(i);
           in_grd_2(i,in_grd_2(i,:)>grd(i))=grd(i);
        end
        for i =1:size(i_pos,2)
            BoxPhantom(in_grd_1(1,i):in_grd_2(1,i),in_grd_1(2,i):in_grd_2(2,i),...
                in_grd_1(3,i):in_grd_2(3,i)) = i_val(i);
        end
    end
end