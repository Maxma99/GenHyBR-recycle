% v1
function  obj = detectorSettingSTIFT_1(obj)
    %step 1: calculate the orientation of the detector
    %step 2: calculate the 2D position of the detector
    %step 3: calculate the 3rd dimension of the detector position
        
    %% step 1: calculate the orientation of the detector
    pos_offset = 1; % an offset value to calculate the 3rd dimension
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
    
   
    
    %%  step 2: calculate the 2D position of the detector
    switch obj.machine_mode % machine mode 
        case {1,3} % machine mode 
            switch detector_orientation % detection orientation
                case 1 % -x: which means the highest position in x axis
                        t_det_number = prod(obj.detector.d_number);
                        det_y_range = [min(obj.node(:,2))+obj.detector.d_edge ...
                            max(obj.node(:,2))-obj.detector.d_edge];
                        det_z_range = [min(obj.node(:,3))+obj.detector.d_edge ...
                            max(obj.node(:,3))-obj.detector.d_edge];
                        [det_y,det_z]=meshgrid...
                            (linspace(det_y_range(1),det_y_range(2),obj.detector.d_number(1)),...
                            linspace(det_z_range(1),det_z_range(2),obj.detector.d_number(2)));
                        list_yz = [reshape(det_y,[],1) reshape(det_z,[],1)];

                case 2 % -y: which means the highest position in y axis
                        t_det_number = prod(obj.detector.d_number);
                        det_x_range = [min(obj.node(:,1))+obj.detector.d_edge ...
                            max(obj.node(:,1))-obj.detector.d_edge];
                        det_z_range = [min(obj.node(:,3))+obj.detector.d_edge ...
                            max(obj.node(:,3))-obj.detector.d_edge];
                        [det_x,det_z]=meshgrid...
                            (linspace(det_x_range(1),det_x_range(2),obj.detector.d_number(1)),...
                            linspace(det_z_range(1),det_z_range(2),obj.detector.d_number(2)));
                        list_xz = [reshape(det_x,[],1) reshape(det_z,[],1)];

                case 3 % -z: which means the highest position in z axis
                        t_det_number = prod(obj.detector.d_number);
                        det_x_range = [min(obj.node(:,1))+obj.detector.d_edge ...
                            max(obj.node(:,1))-obj.detector.d_edge];
                        det_y_range = [min(obj.node(:,2))+obj.detector.d_edge ...
                            max(obj.node(:,2))-obj.detector.d_edge];
                        [det_x,det_y]=meshgrid...
                            (linspace(det_x_range(1),det_x_range(2),obj.detector.d_number(1)),...
                            linspace(det_y_range(1),det_y_range(2),obj.detector.d_number(2)));
                        list_xy = [reshape(det_x,[],1) reshape(det_y,[],1)];

                case 4 % +x: which means the lowest position in x axis
                        t_det_number = prod(obj.detector.d_number);
                        det_y_range = [min(obj.node(:,2))+obj.detector.d_edge ...
                            max(obj.node(:,2))-obj.detector.d_edge];
                        det_z_range = [min(obj.node(:,3))+obj.detector.d_edge ...
                            max(obj.node(:,3))-obj.detector.d_edge];
                        [det_y,det_z]=meshgrid...
                            (linspace(det_y_range(1),det_y_range(2),obj.detector.d_number(1)),...
                            linspace(det_z_range(1),det_z_range(2),obj.detector.d_number(2)));
                        list_yz = [reshape(det_y,[],1) reshape(det_z,[],1)];

                case 5 % +y: which means the lowest position in y axis
                        t_det_number = prod(obj.detector.d_number);
                        det_x_range = [min(obj.node(:,1))+obj.detector.d_edge ...
                            max(obj.node(:,1))-obj.detector.d_edge];
                        det_z_range = [min(obj.node(:,3))+obj.detector.d_edge ...
                            max(obj.node(:,3))-obj.detector.d_edge];
                        [det_x,det_z]=meshgrid...
                            (linspace(det_x_range(1),det_x_range(2),obj.detector.d_number(1)),...
                            linspace(det_z_range(1),det_z_range(2),obj.detector.d_number(2)));
                        list_xz = [reshape(det_x,[],1) reshape(det_z,[],1)];

                case 6 % +z: which means the lowest position in z axis
                        t_det_number = prod(obj.detector.d_number);
                        det_x_range = [min(obj.node(:,1))+obj.detector.d_edge ...
                            max(obj.node(:,1))-obj.detector.d_edge];
                        det_y_range = [min(obj.node(:,2))+obj.detector.d_edge ...
                            max(obj.node(:,2))-obj.detector.d_edge];
                        [det_x,det_y]=meshgrid...
                            (linspace(det_x_range(1),det_x_range(2),obj.detector.d_number(1)),...
                            linspace(det_y_range(1),det_y_range(2),obj.detector.d_number(2)));
                        list_xy = [reshape(det_x,[],1) reshape(det_y,[],1)];
            end % detection orientation
        case 2 % machine mode
            switch obj.detector.d_mode % detector setting mode 
                case 1 % detector setting mode 
                case 2 % detector setting mode 
            end % detector setting mode 
    end % machine mode 
    %%  step 3: calculate the 3rd dimension of the detector position
    switch obj.machine_mode % machine mode 
        case {1,3} % machine mode
            switch obj.detector.d_mode % detector setting mode 
                case 1 % detector setting mode
                    switch detector_orientation % detection orientation
                        case 1 % -x: which means the highest position in x axis
                                list_x = zeros(t_det_number,1);
                                x0 = max(obj.node(:,1));
                                %%%construction%%%%%%%%%%%%%%%%%%%%%%%%
                                for i = 1: t_det_number
                                    list_x(i)=x0;
                                end
                                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                d_virtualPos = [list_x list_yz];

                        case 2 % -y: which means the highest position in y axis
                                list_y = zeros(t_det_number,1);
                                y0 = max(obj.node(:,2));
                                %%%construction%%%%%%%%%%%%%%%%%%%%%%%%
                                for i = 1: t_det_number
                                        list_y(i)=y0;
                                end
                                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                d_virtualPos =[list_xz(:,1) list_y list_xz(:,2)];                        

                        case 3 % -z: which means the highest position in z axis
                                list_z = zeros(t_det_number,1);
                                z0 = max(obj.node(:,3));
                                %%%construction%%%%%%%%%%%%%%%%%%%%%%%%
                                for i = 1: t_det_number
                                        list_z(i)=z0;
                                end
                                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                d_virtualPos =[list_xy(:,1) list_xy(:,2) list_z]; 

                        case 4 % +x: which means the lowest position in x axis
                                list_x = zeros(t_det_number,1);
                                x0 = min(obj.node(:,1));
                                %%%construction%%%%%%%%%%%%%%%%%%%%%%%%
                                for i = 1: t_det_number
                                    list_x(i)=x0;
                                end
                                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                d_virtualPos = [list_x list_yz];

                        case 5 % +y: which means the lowest position in y axis
                                list_y = zeros(t_det_number,1);
                                y0 = min(obj.node(:,2));
                                %%%construction%%%%%%%%%%%%%%%%%%%%%%%%
                                for i = 1: t_det_number
                                       list_y(i)=y0;
                                end
                                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                d_virtualPos =[list_xz(:,1) list_y list_xz(:,2)];

                        case 6 % +z: which means the lowest position in z axis
                                list_z = zeros(t_det_number,1);
                                z0 = min(obj.node(:,3));
                                %%%construction%%%%%%%%%%%%%%%%%%%%%%%%
                                for i = 1: t_det_number
                                        list_z(i)=z0;
                                end
                                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                d_virtualPos =[list_xy(:,1) list_xy(:,2) list_z];
                    end % detection orientation
                case 2 % detector setting mode 
                    switch  detector_orientation % detection orientation
                        case 1 % -x: which means the highest position in x axis
                                list_x = zeros(t_det_number,1);
                                x0 = max(obj.node(:,1))+pos_offset;
                                %%%construction%%%%%%%%%%%%%%%%%%%%%%%%
                                for i = 1: t_det_number
                                    p0=[x0 list_yz(i,1) list_yz(i,2)];
                                    [~,~,~,~,xnode]=raysurf(p0,obj.detector.d_normal_dA,obj.node,obj.face);
                                    temp_node=max(xnode(:,1));
                                    if temp_node>0
                                        list_x(i) = temp_node;
                                    else
                                        list_x(i) = x0;
                                    end
                                end
                                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                d_virtualPos = [list_x list_yz];

                        case 2 % -y: which means the highest position in y axis
                                list_y = zeros(t_det_number,1);
                                y0 = max(obj.node(:,2))+pos_offset;
                                %%%construction%%%%%%%%%%%%%%%%%%%%%%%%
                                for i = 1: t_det_number
                                        p0=[list_xz(i,1) y0 list_xz(i,2)];
                                        [~,~,~,~,ynode]=raysurf(p0,obj.detector.d_normal_dA,obj.node,obj.face);
                                        temp_node=max(ynode(:,2));
                                        if temp_node>0
                                            list_y(i) = temp_node;
                                        else
                                            list_y(i) = y0;
                                        end
                                end
                                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                d_virtualPos =[list_xz(:,1) list_y list_xz(:,2)];                          

                        case 3 % -z: which means the highest position in z axis
                                list_z = zeros(t_det_number,1);
                                z0 = max(obj.node(:,3))+pos_offset;
                                %%%construction%%%%%%%%%%%%%%%%%%%%%%%%
                                for i = 1: t_det_number
                                        p0=[list_xy(i,1) list_xy(i,2) z0 ];
                                        [~,~,~,~,znode]=raysurf(p0,obj.detector.d_normal_dA,obj.node,obj.face);
                                        temp_node=max(znode(:,3));
                                        if temp_node>0
                                            list_z(i) = temp_node;
                                        else
                                            list_z(i) = z0;
                                        end
                                end
                                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                d_virtualPos =[list_xy(:,1) list_xy(:,2) list_z];

                        case 4 % +x: which means the lowest position in x axis
                                list_x = zeros(t_det_number,1);
                                x0 = min(obj.node(:,1))-pos_offset;
                                %%%construction%%%%%%%%%%%%%%%%%%%%%%%%
                                for i = 1: t_det_number
                                    p0=[x0 list_yz(i,1) list_yz(i,2)];
                                    [~,~,~,~,xnode]=raysurf(p0,obj.detector.d_normal_dA,obj.node,obj.face);
                                    temp_node=min(xnode(:,1));
                                    if temp_node>0
                                        list_x(i) = temp_node;
                                    else
                                        list_x(i) = x0;
                                    end
                                end
                                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                d_virtualPos = [list_x list_yz];

                        case 5 % +y: which means the lowest position in y axis
                                list_y = zeros(t_det_number,1);
                                y0 = min(obj.node(:,2))-pos_offset;
                                %%%construction%%%%%%%%%%%%%%%%%%%%%%%%
                                for i = 1: t_det_number
                                        p0=[list_xz(i,1) y0 list_xz(i,2)];
                                        [~,~,~,~,ynode]=raysurf(p0,obj.detector.d_normal_dA,obj.node,obj.face);
                                        temp_node=min(ynode(:,2));
                                        if temp_node>0
                                            list_y(i) = temp_node;
                                        else
                                            list_y(i) = y0;
                                        end
                                end
                                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                d_virtualPos =[list_xz(:,1) list_y list_xz(:,2)]; 

                        case 6 % +z: which means the lowest position in z axis
                                list_z = zeros(t_det_number,1);
                                z0 = min(obj.node(:,3))-pos_offset;
                                %%%construction%%%%%%%%%%%%%%%%%%%%%%%%
                                for i = 1: t_det_number
                                        p0=[list_xy(i,1) list_xy(i,2) z0 ];
                                        [~,~,~,~,znode]=raysurf(p0,obj.detector.d_normal_dA,obj.node,obj.face);
                                        temp_node=min(znode(:,3));
                                        if temp_node>0
                                            list_z(i) = temp_node;
                                        else
                                            list_z(i) = z0;
                                        end
                                end
                                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                d_virtualPos =[list_xy(:,1) list_xy(:,2) list_z];
                    end  % detection orientation
            end % detector setting mode 
        case 2 % machine mode
            switch obj.detector.d_mode % detector setting mode 
                case 1 % detector setting mode 
                case 2 % detector setting mode 
            end % detector setting mode 
    end % machine mode  
    %% mask circle DOT-FMT phantom
%     center = (obj.detector.d_number+1)/2;
%     rad = center(1)-1;
%     mask = zeros(obj.detector.d_number);
%     for i = 1:obj.detector.d_number(1)
%         for j = 1:obj.detector.d_number(2)
%             if(i-center(1))^2 + (j-center(2))^2 <= rad^2
%                 mask(i,j) = 1;
%             end
%         end
%     end
%     d_virtualPos_temp = zeros(sum(mask(:)),3);
%     cnt = 0;
%     for i = 1:size(d_virtualPos,1)
%         if(mask(ceil(d_virtualPos(i,1))+1,ceil(d_virtualPos(i,2))+1))
%             cnt = cnt + 1;
%             d_virtualPos_temp(cnt,:) = d_virtualPos(i,:);
%         end
%     end
%     d_virtualPos = d_virtualPos_temp;

    %%
    obj.detector.d_virtualPos = d_virtualPos;
    
    
end