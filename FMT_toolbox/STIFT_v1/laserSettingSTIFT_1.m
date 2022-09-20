% v1
function obj = laserSettingSTIFT_1(obj)

    thre = 0.4; % threshold to defining laser points
    %% set the laser orientation 光照方向 从左往右 从上往下 从前往后等6种
    if (obj.laser.l_direction(1) <0)... 
            &&(obj.laser.l_direction(2)==0)...
            &&(obj.laser.l_direction(3)==0)
        laser_orientation = 1 ;
    end
    if (obj.laser.l_direction(2) <0)...
            &&(obj.laser.l_direction(1)==0)...
            &&(obj.laser.l_direction(3)==0)
        laser_orientation = 2 ;
    end
    if (obj.laser.l_direction(3) <0)...
            &&(obj.laser.l_direction(1)==0)...
            &&(obj.laser.l_direction(2)==0)
        laser_orientation = 3 ;
    end
    if (obj.laser.l_direction(1) >0)...
            &&(obj.laser.l_direction(2)==0)...
            &&(obj.laser.l_direction(3)==0)
        laser_orientation = 4 ;
    end
    if (obj.laser.l_direction(2) >0)...
            &&(obj.laser.l_direction(1)==0)...
            &&(obj.laser.l_direction(3)==0)
        laser_orientation = 5 ;
    end
    if (obj.laser.l_direction(3) >0)...
            &&(obj.laser.l_direction(1)==0)...
            &&(obj.laser.l_direction(2)==0)
        laser_orientation = 6 ;
    end
    
    switch laser_orientation
        case 1 % -x;
            if  strcmp(obj.laser.l_design,'regu')
                t_laser_number = prod(obj.laser.l_number);
                scan_y_range = [obj.laser.l_region(1,1) obj.laser.l_region(1,2)];
                scan_z_range = [obj.laser.l_region(2,1) obj.laser.l_region(2,2)];
                [source_y,source_z]=meshgrid...
                    (linspace(scan_y_range(1),scan_y_range(2),obj.laser.l_number(1)),...
                    linspace(scan_z_range(1),scan_z_range(2),obj.laser.l_number(2)));
                list_yz = [reshape(source_y,[],1) reshape(source_z,[],1)];
                list_x = zeros(t_laser_number,1);
            elseif strcmp(obj.laser.l_design,'load')
                %[list_yz t_laser_number Ex_handle] = getLaserPosition(obj.laser.l_pos_used);
                %list_yz = [list_yz(:,1)-obj.measurement.origin(1) ...
                %    list_yz(:,2)-obj.measurement.origin(2)]*obj.measurement.pixel_size;
                %list_x = zeros(t_laser_number,1);  
                disp('Currently loading is only possible for +/-z beam');   
            else
                disp('wrong laser settings!');    
            end
            x0 = max(obj.node(:,1))+1;
            %%%construction
            for i = 1: t_laser_number
                p0=[x0 list_yz(i,1) list_yz(i,2)];
                [~,~,~,~,xnode]=raysurf(p0,obj.laser.l_direction,obj.node,obj.face);
                list_x(i)=max(xnode(:,1));
            end
            %%%%
            obj.laser.l_position = [list_x list_yz];
            l_index_A = obj.laser.l_position(:,1)<=x0;
            temp_pos = obj.laser.l_position(l_index_A,:);
            clear obj.laser.l_position
            obj.laser.l_position = temp_pos;
            obj.laser.l_pos_used(~l_index_A )=[]; 
        case 2 % -y
            if strcmp(obj.laser.l_design,'regu')
                t_laser_number = prod(obj.laser.l_number);
                scan_x_range = [obj.laser.l_region(1,1) obj.laser.l_region(1,2)];
                scan_z_range = [obj.laser.l_region(2,1) obj.laser.l_region(2,2)];
                [source_x,source_z]=meshgrid...
                    (linspace(scan_x_range(1),scan_x_range(2),obj.laser.l_number(1)),...
                    linspace(scan_z_range(1),scan_z_range(2),obj.laser.l_number(2)));
                list_xz = [reshape(source_x,[],1) reshape(source_z,[],1)];
                list_y = zeros(t_laser_number,1);
            elseif strcmp(obj.laser.l_design,'load')
                %[list_xz t_laser_number Ex_handle] = getLaserPosition(obj.laser.l_pos_used);
                %list_xz = [list_xz(:,1)-obj.measurement.origin(1) ...
                %    list_xz(:,2)-obj.measurement.origin(2)]*obj.measurement.pixel_size;
                %list_y = zeros(t_laser_number,1);
                disp('Currently loading is only possible for +/-z beam');   
            else
                disp('wrong laser settings!');    
            end
            y0=max(obj.node(:,2))+1;
            %%%construction
            for i = 1: t_laser_number
                p0=[list_xz(i,1) y0 list_xz(i,2)];
                [~,~,~,~,ynode]=raysurf(p0,obj.laser.l_direction,obj.node,obj.face);
                list_y(i)=max(ynode(:,2));
            end
            %%%%
            obj.laser.l_position = [list_xz(:,1) list_y list_xz(:,2)];
            l_index_A = obj.laser.l_position(:,2)<=y0;
            temp_pos = obj.laser.l_position(l_index_A,:);
            clear obj.laser.l_position
            obj.laser.l_position = temp_pos;      
            obj.laser.l_pos_used(~l_index_A )=[]; 
        case 3 % -z
            if  strcmp(obj.laser.l_design,'regu')
                t_laser_number = prod(obj.laser.l_number);
                scan_x_range = [obj.laser.l_region(1,1) obj.laser.l_region(1,2)];
                scan_y_range = [obj.laser.l_region(2,1) obj.laser.l_region(2,2)];
                [source_x,source_y]=meshgrid...
                    (linspace(scan_x_range(1),scan_x_range(2),obj.laser.l_number(1)),...
                    linspace(scan_y_range(1),scan_y_range(2),obj.laser.l_number(2)));
            list_xy = [reshape(source_x,[],1) reshape(source_y,[],1)];
            list_z = zeros(t_laser_number,1);
            elseif strcmp(obj.laser.l_design,'load')
                [list_xy,t_laser_number] = loadLaser(thre,obj.project_directory,obj.laser.l_pos_used);
                list_xy =[list_xy(:,2)-obj.measurement.origin(1) ...
                	list_xy(:,1)-obj.measurement.origin(2)].*obj.measurement.pixel_size;
                list_z = zeros(t_laser_number,1);                  
            else
                disp('wrong laser settings!');    
            end
            z0=max(obj.node(:,3))+1;
            %%%construction
            for i = 1: t_laser_number
                p0=[list_xy(i,1) list_xy(i,2) z0 ];
                [~,~,~,~,znode]=raysurf(p0,obj.laser.l_direction,obj.node,obj.face);
                list_z(i)=max(znode(:,3));
            end
            %%%%
            obj.laser.l_position = [list_xy list_z];
            l_index_A = obj.laser.l_position(:,3)<=z0;
            temp_pos = obj.laser.l_position(l_index_A,:);
            clear obj.laser.l_position
            obj.laser.l_position = temp_pos;   
            obj.laser.l_pos_used(~l_index_A )=[]; 
        case 4 % +x;
            if  strcmp(obj.laser.l_design,'regu')
                t_laser_number = prod(obj.laser.l_number);
                scan_y_range = [obj.laser.l_region(1,1) obj.laser.l_region(1,2)];
                scan_z_range = [obj.laser.l_region(2,1) obj.laser.l_region(2,2)];
                [source_y source_z]=meshgrid...
                    (linspace(scan_y_range(1),scan_y_range(2),obj.laser.l_number(1)),...
                    linspace(scan_z_range(1),scan_z_range(2),obj.laser.l_number(2)));
                list_yz = [reshape(source_y,[],1) reshape(source_z,[],1)];
                list_x = zeros(t_laser_number,1);
            elseif strcmp(obj.laser.l_design,'load')
                %[list_yz t_laser_number Ex_handle] = getLaserPosition(obj.laser.l_pos_used);
                %list_yz = [list_yz(:,1)-obj.measurement.origin(1) ...
                %    list_yz(:,2)-obj.measurement.origin(2)]*obj.measurement.pixel_size;
                %list_x = zeros(t_laser_number,1); 
                disp('Currently loading is only possible for +/-z beam');   
            else
                disp('wrong laser settings!');    
            end
            x0 = min(obj.node(:,1))-1;
            %%%construction
            for i = 1: t_laser_number
                p0=[x0 list_yz(i,1) list_yz(i,2)];
                [~,~,~,~,xnode]=raysurf(p0,obj.laser.l_direction,obj.node,obj.face);
                list_x(i)=min(xnode(:,1));
            end
            %%%%
            obj.laser.l_position = [list_x list_yz];
            l_index_A = obj.laser.l_position(:,1)>=x0;
            temp_pos = obj.laser.l_position(l_index_A,:);
            clear obj.laser.l_position
            obj.laser.l_position = temp_pos;
            obj.laser.l_pos_used(~l_index_A )=[]; 
        case 5% +y
            if  strcmp(obj.laser.l_design,'regu')
                t_laser_number = prod(obj.laser.l_number);
                scan_x_range = [obj.laser.l_region(1,1) obj.laser.l_region(1,2)];
                scan_z_range = [obj.laser.l_region(2,1) obj.laser.l_region(2,2)];
                [source_x,source_z]=meshgrid...
                    (linspace(scan_x_range(1),scan_x_range(2),obj.laser.l_number(1)),...
                    linspace(scan_z_range(1),scan_z_range(2),obj.laser.l_number(2)));
                list_xz = [reshape(source_x,[],1) reshape(source_z,[],1)];
                list_y = zeros(t_laser_number,1);
            elseif strcmp(obj.laser.l_design,'load')
                %[list_xz t_laser_number Ex_handle] = getLaserPosition(obj.laser.l_pos_used);
                %list_xz = [list_xz(:,1)-obj.measurement.origin(1) ...
                %    list_xz(:,2)-obj.measurement.origin(2)]*obj.measurement.pixel_size;
                %list_y = zeros(t_laser_number,1);  
                disp('Currently loading is only possible for +/-z beam');   
            else
                disp('wrong laser settings!');    
            end
            y0=min(obj.node(:,2))-1;
            %%%construction
            for i = 1: t_laser_number
                p0=[list_xz(i,1) y0 list_xz(i,2)];
               [~,~,~,~,ynode]=raysurf(p0,obj.laser.l_direction,obj.node,obj.face);
                list_y(i)=min(ynode(:,2));
            end
            %%%%
            obj.laser.l_position = [list_xz(:,1) list_y list_xz(:,2)];
            l_index_A = obj.laser.l_position(:,2)>=y0;
            temp_pos = obj.laser.l_position(l_index_A,:);
            clear obj.laser.l_position
            obj.laser.l_position = temp_pos;     
            obj.laser.l_pos_used(~l_index_A )=[]; 
        case 6  % +z
            if  strcmp(obj.laser.l_design,'regu')
                t_laser_number = prod(obj.laser.l_number);
                scan_x_range = [obj.laser.l_region(1,1) obj.laser.l_region(1,2)];
                scan_y_range = [obj.laser.l_region(2,1) obj.laser.l_region(2,2)];
                [source_x,source_y]=meshgrid...
                    (linspace(scan_x_range(1),scan_x_range(2),obj.laser.l_number(1)),...
                    linspace(scan_y_range(1),scan_y_range(2),obj.laser.l_number(2)));
                list_xy = [reshape(source_x,[],1) reshape(source_y,[],1)];
                list_z = zeros(t_laser_number,1);
            elseif strcmp(obj.laser.l_design,'load')
                %[list_xy t_laser_number Ex_handle] = getLaserPosition(obj.laser.l_pos_used);
                %list_xy = [list_xy(:,1)-obj.measurement.origin(1) ...
                %    list_xy(:,2)-obj.measurement.origin(2)]*obj.measurement.pixel_size;
                %list_z = zeros(t_laser_number,1);   
                [list_xy,t_laser_number] = loadLaser(thre,obj.project_directory,obj.laser.l_pos_used);
                list_xy =[list_xy(:,2)-obj.measurement.origin(1) ...
                	list_xy(:,1)-obj.measurement.origin(2)].*obj.measurement.pixel_size;
                list_z = zeros(t_laser_number,1); 
            else
                disp('wrong laser settings!');    
            end
            z0=min(obj.node(:,3))-1;
            %%%construction
            for i = 1: t_laser_number
                p0=[list_xy(i,1) list_xy(i,2) z0 ];
               [~,~,~,~,znode]=raysurf(p0,obj.laser.l_direction,obj.node,obj.face);
                list_z(i)=min(znode(:,3));
            end
            %%%%
            obj.laser.l_position = [list_xy list_z];
            l_index_A = obj.laser.l_position(:,3)>=z0;
            temp_pos = obj.laser.l_position(l_index_A,:);
            clear obj.laser.l_position
           obj.laser.l_position = temp_pos;     
%            obj.laser.l_position = load("l_position.mat");
            obj.laser.l_pos_used(~l_index_A )=[]; 
    end
    x_range = [0.15,0.85];
y_range = [0.15,0.85];
    nx = 10;
    ny = 10;
    [det_x,det_z] = meshgrid(linspace(x_range(1),x_range(2),nx),linspace(y_range(1),y_range(2),ny));
    scan_points_xy = [reshape(det_x,[],1) reshape(det_z,[],1)];
    scan_points_xy_1 = [scan_points_xy,0.01*ones(size(scan_points_xy,1),1)];
    obj.laser.l_position = obj.phantom_parameter.x * scan_points_xy_1;
end