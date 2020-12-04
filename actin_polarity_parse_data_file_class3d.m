function [alg_struct] = actin_polarity_parse_data_file_class3d(relion_data_file)
%%%%%%%%%
%%%%%%%%%
%

% data_
% 
% loop_ 
% _rlnImageName #1 
% _rlnMicrographName #2 
% _rlnAngleRotPrior #3 
% _rlnAngleTiltPrior #4 
% _rlnAnglePsiPrior #5 
% _rlnOriginXPrior #6 
% _rlnOriginYPrior #7 
% _rlnGroupNumber #8 
% _rlnAngleRot #9 
% _rlnAngleTilt #10 
% _rlnAnglePsi #11 
% _rlnOriginX #12 
% _rlnOriginY #13 
% _rlnClassNumber #14 
% _rlnNormCorrection #15 
% _rlnRandomSubset #16 
% _rlnLogLikeliContribution #17 
% _rlnMaxValueProbDistribution #18 
% _rlnNrOfSignificantSamples #19 

% Open particles list
fid = fopen(relion_data_file);

% Parse particles list
zaehler = 1;
while 1
    
    tline = fgetl(fid);
    
    values = regexp(tline,' *','split');
    
    if tline==-1 | size(values,2) < 19
         disp('End of file reached!');
         break;
    end
    
    alg_struct(zaehler).data_line = values;
    
    alg_struct(zaehler).image_name = values{1};
    image_name_pre = alg_struct(zaehler).image_name;
    alg_struct(zaehler).particle_indx = str2num(image_name_pre(1:6));
    alg_struct(zaehler).micrograph_name = values{2};
    micrograph_name_pre = alg_struct(zaehler).micrograph_name;
    alg_struct(zaehler).micrograph_indx = str2num(micrograph_name_pre(21:21));
    
    alg_struct(zaehler).angle_rot = str2double(values{9});
    alg_struct(zaehler).angle_tilt = str2double(values{10});
    alg_struct(zaehler).angle_psi = str2double(values{11});
    alg_struct(zaehler).origin_x = str2double(values{12});
    alg_struct(zaehler).origin_y = str2double(values{13});
    
    alg_struct(zaehler).class_indx = str2double(values{14});
    
    zaehler = zaehler + 1;
    
end

% Close particles list
fclose(fid);

