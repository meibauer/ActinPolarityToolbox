function [image_name,particle_indx,angle_rot,angle_tilt,angle_psi,origin_x,origin_y,class_indx,rlnLogLikeliContribution,rlnMaxValueProbDistribution,rlnNrOfSignificantSamples] = actin_polarity_parse_data_file_refine3d(relion_3d_data_file_mod)
%%%%%%%%
%%%%%%%%%
%
% This function reads a modified Relion Refine3D data file.
%
% INPUT
% relion_3d_data_file_mod --- Name of the modified Relion Refine3D data file 
%
% OUTPUT
% particle_list --- List of particle images
% particle_indx --- Vector of particle indices
% angle_rot --- Relion rot angle
% angle_tilt --- Relion tilt angle
% angle_psi --- Relion psi angle
% origin_x --- Vector of in-plane translations in x-direction
% origin_y --- Vector of in-plane translations in y-direction

% Open particles list
fid = fopen(relion_3d_data_file_mod);

% Parse particles list
zaehler = 1;
while 1
    
    tline = fgetl(fid);
    
    if(tline==-1)
         break;
    end
    
    if regexp(tline,'.mrc')
    
    values = regexp(tline,' *','split');
    
    image_name{zaehler} = values{1};
    
    angle_rot(zaehler) = str2double(values{9});    % angle_rot
    angle_tilt(zaehler) = str2double(values{10});   % angle_tilt
    angle_psi(zaehler) = str2double(values{11});    % angle_psi
    
    origin_x(zaehler) = str2double(values{12});
    origin_y(zaehler) = str2double(values{13});
    
    class_indx(zaehler) = str2double(values{14});
    
    rlnLogLikeliContribution(zaehler) = str2double(values{17});
    rlnMaxValueProbDistribution(zaehler) = str2double(values{18});
    rlnNrOfSignificantSamples(zaehler) = str2double(values{19});
    
    zaehler = zaehler + 1;
    
    end

end

% Close particles list
fclose(fid);

% Extract particle indx
for k=1:size(image_name,2)
    image_name_line = image_name{k};
    particle_indx(k) = str2num(image_name_line(1:6));
end

