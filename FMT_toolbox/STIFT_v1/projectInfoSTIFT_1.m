% v1
function  obj  = projectInfoSTIFT_1(obj)
        default_folder_name = pwd;
        fprintf(1,'Create a FMT project ...\n');
        obj.project_directory = uigetdir(default_folder_name,...
            'Please select a project folder');% folder to save the project
        fprintf(1,'Choose the data buffer folder ...\n');
        obj.data_buffer_directory = uigetdir(default_folder_name,...
            'Please data buffer folder');% folder to save the data buffer
        obj.date = 20220221; % experiment date
        obj.time = 0; % experiment time
        obj.object_ID = 'FMT test'; % the object ID (mouse ID)
        obj.user_name = 'FMT guy'; % the name of user
        if obj.machine_mode == 3
            %% read from the standalone FMT
            % obj.measurement.pixel_size = 60/222;% typical value for
            % 1)standalone FMT (mm): 60/222; 2) MRI:20/128 
            % obj.measurement.origin = [x0 y0];
            obj.measurement.pixel_size = 54/660 ; %mm
            obj.measurement.origin = [0 0];
            
%             fprintf(1,'loading an excitation map ...\n');
%             [obj.measurement.ex_file,obj.measurement.ex_folder] = uigetfile('*.tif','Please select an excitation map',...
%             obj.project_directory);
%             if isequal(obj.measurement.ex_file,0)
%                 disp('no excitation selected')
%             else
%                 disp('done')
%             end
% 
%             fprintf(1,'loading an emission map ...\n');
%             [obj.measurement.em_file,obj.measurement.em_folder] = uigetfile('*.tif','Please select an emission map',...
%             obj.measurement.ex_folder);
%             if isequal(obj.measurement.em_file,0)
%                 disp('no emission selected')
%             else
%                 disp('done')
%             end
% 
%             fprintf(1,'loading a whitelight image ...\n');
%             [obj.measurement.wt_file,obj.measurement.wt_folder] = uigetfile('*.tif','Please select wt background_',...
%             obj.measurement.ex_folder);
%             if isequal(obj.measurement.wt_file,0)
%                 disp('no background selected')
%             else
%                 disp('done')
%             end
        end
end