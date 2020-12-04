% Actin polarity / averaging / box size 144 pixels 

%% Create particle list and alignment structure

% Generate BasePath
actin_polarity_construct_base_path;

% Define microscope parameter
Objectpixelsize = 0.3443;%nm

% Define particle path / bin0_144
ParticlesPath = '/home/Medalia/scratchX03/ActinPolarity/particles/actin_particles_bin0_144';

% Loop over actin particles
zaehler = 1;
zaehler_particles = 0;

for k=1:size(BasePath,2)
    
    % Check if particles are selected
    try
         load([BasePath{k} '/cor/actin_polarity_cor_final_144.mat'],'plist');
         particles_selected = 1;
         
    catch
         particles_selected = 0;
    end
    
    % Perform task
    if particles_selected
        
         for p=1:size(plist,1)
             
              % Identify particle
              actin_particle_pre = [ParticlesPath '/actin_particle_t_' num2str(k) '_p_' num2str(p) '.em'];
              disp(actin_particle_pre);
              
              % Create new particle identifier
              actin_particle = [ParticlesPath '/actin_particle_' actin_polarity_construct_particle_identifier(k,p,zaehler) '.em'];
              disp(actin_particle);
              
              % Rename particle
              unix(['mv ' actin_particle_pre ' ' actin_particle]);
              
              particles_list_pre{zaehler} = actin_particle_pre;
              particles_list{zaehler} = actin_particle;
              
              %disp(zaehler);
              
              zaehler = zaehler + 1;
             
         end
         
         zaehler_particles = zaehler_particles + size(plist,1);
         
         clear plist;
         
    end
    
    clear particles_selected;
    
    %disp(k);
    
end

% Build particle list
zaehler = 1;
for k=1:size(particles_list,2)
    
    particles_list_combined{zaehler} = particles_list{k};
    zaehler = zaehler + 1;
    
    particles_list_combined{zaehler} = 'singleaxiswedge 30 72';
    zaehler = zaehler + 1;
    
    particles_list_combined{zaehler} = 'allpass';
    zaehler = zaehler + 1;
    
    particles_list_combined{zaehler} = ' ';
    zaehler = zaehler + 1;
    
end

% Write empty file
placeh = num2str(0);
          save('/home/Medalia/Projects7/Bruno/ActinPolarity/averaging/actin_polarity_particles_list_bin0_144.txt', 'placeh', '-ASCII');

% Write particles into file
fid = fopen('/home/Medalia/Projects7/Bruno/ActinPolarity/averaging/actin_polarity_particles_list_bin0_144.txt', 'w');
for k=1:size(particles_list_combined,2)
    fprintf(fid, '%s\n', char(particles_list_combined{k}));        
end
fclose(fid);

% Save workspace
save('/home/Medalia/Projects7/Bruno/ActinPolarity/averaging/workspace_actin_polarity_averaging_144_step_1.mat');
clear all;


%% Project particles

% Parse particles list
[particles_list] = ParseParticleList('/home/Medalia/Projects7/Bruno/ActinPolarity/averaging/actin_polarity_particles_list_bin0_144.txt');

% Define objectpixelsize
Objectpixelsize = 0.3443 .* 10;%A

% Calculate resolution
f = (1:72)./(2.*Objectpixelsize.*72);
x = round((1./f).*10)./10;%A

% Project particles
for k=1:size(particles_list,2)
    
    % Load particle
    particle = double(tom_emreadc3f(particles_list{k}));
    
    % Project actin filament
    proj = mean(particle(:,:,72-16+1:72+16),3);%---> Projection thickness = 11 nm (32 pixels)
    
    % Normalize particle
    proj = (proj - mean(proj(:)))./std(proj(:));
    
    % Dissect file name into parts
    [~,particle_name] = fileparts(particles_list{k});
    
    % Write to disk
    tom_mrcwrite(single(proj),'name',['/home/Medalia/scratchX03/ActinPolarity/particles/actin_particles_bin0_144_proj/' particle_name '.mrc']);
    
                      % Prepare proj list
                      proj_list{k} = ['/home/Medalia/scratchX03/ActinPolarity/particles/actin_particles_bin0_144_proj/' particle_name '.mrc' '   ' '/home/Medalia/scratchX03/ActinPolarity/particles/actin_particles_bin0_144_proj/'  particle_name(1:21) '.mrc'];
                      
    disp(k);
    
end

% Write proj list
placeh = num2str(0);
          save('/home/Medalia/Projects7/Bruno/ActinPolarity/averaging/actin_polarity_particles_list_bin0_144_proj.star', 'placeh', '-ASCII');
fid = fopen('/home/Medalia/Projects7/Bruno/ActinPolarity/averaging/actin_polarity_particles_list_bin0_144_proj.star', 'w');
for k=1:size(proj_list,2)
    fprintf(fid, '%s\n', char(proj_list{k}));
end
fclose(fid);

% Add the following header:

% data_
% loop_
% _rlnImageName
% _rlnMicrographName

% Normalize particles in Relion
%mpirun --np 32 relion_preprocess_mpi --operate_on /home/Medalia/scratchX03/ActinPolarity/particles/actin_polarity_particles_list_bin0_144_proj.star --norm --bg_radius 72 --invert_contrast --operate_out /home/Medalia/scratchX03/ActinPolarity/particles/actin_polarity_particles_list_bin0_144_proj_norm_inv --set_angpix 3.443

% Save workspace
save('/home/Medalia/Projects7/Bruno/ActinPolarity/averaging/workspace_actin_polarity_averaging_144_step_2.mat');
clear all;


%% Prepare 2D classification with priors / Class2D / job001 / it010 / prealignment with template library

% Parse data file / "this is the prealignment"
[alg_struct] = actin_polarity_parse_data_file_class2d_prealg('/home/Medalia/Projects7/Bruno/ActinPolarity/averaging/rln_144/Class2D/job001/run_it010_data_mod.star');

% Create new starfile
zaehler = 1;
for k=1:size(alg_struct,2)
    
    relion_data_file_out{zaehler} = [alg_struct(k).image_name '   ' alg_struct(k).micrograph_name '   ' num2str(0) '   ' num2str(0) '   ' num2str(alg_struct(k).angle_psi) '   ' num2str(alg_struct(k).origin_x) '   ' num2str(alg_struct(k).origin_y)];
    
    zaehler = zaehler + 1;
    
end

% Write starfile
placeh = num2str(0);
       save('/home/Medalia/Projects7/Bruno/ActinPolarity/averaging/rln_144/Class2D/job001/particles_job001_it010_prior_mod.star', 'placeh', '-ASCII');
fid = fopen('/home/Medalia/Projects7/Bruno/ActinPolarity/averaging/rln_144/Class2D/job001/particles_job001_it010_prior_mod.star', 'w');
for k=1:size(relion_data_file_out,2)
    fprintf(fid, '%s\n', char(relion_data_file_out{k}));
end
fclose(fid);

% Add header manually

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


%% Prepare 3D reconstruction with priors / Class2D / job003 / it100 / Select / job004 / alignment with prior psi angle

% Parse data file / "this is the psi finealignment"
[alg_struct] = actin_polarity_parse_data_file_class2d_finealg('/home/Medalia/Projects7/Bruno/ActinPolarity/averaging/rln_144/Select/job004/particles_mod.star');

% Create new starfile
zaehler = 1;
group_indx = 1;
for k=1:size(alg_struct,2)
    
    if group_indx == 1
         relion_data_file_out{zaehler} = [alg_struct(k).image_name '   ' 'actin_polarity_tomo_' num2str(group_indx) '   ' num2str((rand(1)-0.5).*360) '   ' num2str(90) '   ' num2str(alg_struct(k).angle_psi) '   ' num2str(alg_struct(k).origin_x) '   ' num2str(alg_struct(k).origin_y)];
         group_indx = 2;
    else
         relion_data_file_out{zaehler} = [alg_struct(k).image_name '   ' 'actin_polarity_tomo_' num2str(group_indx) '   ' num2str((rand(1)-0.5).*360) '   ' num2str(90) '   ' num2str(alg_struct(k).angle_psi) '   ' num2str(alg_struct(k).origin_x) '   ' num2str(alg_struct(k).origin_y)];
         group_indx = 1;
    end
    
    zaehler = zaehler + 1;
    
end

% Write starfile and add header manually
placeh = num2str(0);
       save('/home/Medalia/Projects7/Bruno/ActinPolarity/averaging/rln_144/Select/job004/particles_job004_prior_mod.star', 'placeh', '-ASCII');
fid = fopen('/home/Medalia/Projects7/Bruno/ActinPolarity/averaging/rln_144/Select/job004/particles_job004_prior_mod.star', 'w');
for k=1:size(relion_data_file_out,2)
    fprintf(fid, '%s\n', char(relion_data_file_out{k}));
end
fclose(fid);

% Add header manually

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


%% Prepare local alignment and classification based on Refine3D / job006

% Parse data file
[alg_struct] = actin_polarity_parse_data_file_class3d('/home/Medalia/Projects7/Bruno/ActinPolarity/averaging/rln_144/Refine3D/job006/run_data_mod.star');

% Create new starfile
zaehler = 1;
group_indx = 1;
for k=1:size(alg_struct,2)
    
    if group_indx == 1
         relion_data_file_out{zaehler} = [alg_struct(k).image_name '   ' 'actin_polarity_tomo_' num2str(group_indx) '   ' num2str(alg_struct(k).angle_rot) '   ' num2str(alg_struct(k).angle_tilt) '   ' num2str(alg_struct(k).angle_psi) '   ' num2str(alg_struct(k).origin_x) '   ' num2str(alg_struct(k).origin_y)];
         group_indx = 2;
    else
         relion_data_file_out{zaehler} = [alg_struct(k).image_name '   ' 'actin_polarity_tomo_' num2str(group_indx) '   ' num2str(alg_struct(k).angle_rot) '   ' num2str(alg_struct(k).angle_tilt) '   ' num2str(alg_struct(k).angle_psi) '   ' num2str(alg_struct(k).origin_x) '   ' num2str(alg_struct(k).origin_y)];
         group_indx = 1;
    end
    
    zaehler = zaehler + 1;
    
end

% Write starfile and add header manually
placeh = num2str(0);
       save('/home/Medalia/Projects7/Bruno/ActinPolarity/averaging/rln_144/Refine3D/job006/particles_job006_prior_mod.star', 'placeh', '-ASCII');
fid = fopen('/home/Medalia/Projects7/Bruno/ActinPolarity/averaging/rln_144/Refine3D/job006/particles_job006_prior_mod.star', 'w');
for k=1:size(relion_data_file_out,2)
    fprintf(fid, '%s\n', char(relion_data_file_out{k}));
end
fclose(fid);

% Add header manually

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


%% Reduce data file / particles_job004_prior_particles_del_translational_priors.star / 20'585 particles ---> 9'931 particles

% Parse data file
[alg_struct] = actin_polarity_parse_data_file_tmp('/home/Medalia/Projects7/Bruno/ActinPolarity/averaging/rln_144/Select/job004/particles_job004_prior_particles_del_translational_priors_mod.star');

% Select 9'931 random indices
indx_sel = randperm(size(alg_struct,2));
indx_sel = sort(indx_sel(1:9931));

% Create new starfile
zaehler = 1;
for k=indx_sel
    
    relion_data_file_out{zaehler} = [alg_struct(k).image_name '   ' alg_struct(k).micrograph_name '   ' num2str(round(alg_struct(k).angle_rot)) '   ' num2str(alg_struct(k).angle_tilt) '   ' num2str(alg_struct(k).angle_psi)];
    
    zaehler = zaehler + 1;
    
end

% Write starfile and add header manually
placeh = num2str(0);
       save('/home/Medalia/Projects7/Bruno/ActinPolarity/averaging/rln_144/Select/job004/particles_job004_prior_particles_del_translational_priors_rand_red_9931_mod.star', 'placeh', '-ASCII');
fid = fopen('/home/Medalia/Projects7/Bruno/ActinPolarity/averaging/rln_144/Select/job004/particles_job004_prior_particles_del_translational_priors_rand_red_9931_mod.star', 'w');
for k=1:size(relion_data_file_out,2)
    fprintf(fid, '%s\n', char(relion_data_file_out{k}));
end
fclose(fid);

% Add header manually

% data_
% 
% loop_ 
% _rlnImageName #1 
% _rlnMicrographName #2 
% _rlnAngleRotPrior #3 
% _rlnAngleTiltPrior #4 
% _rlnAnglePsiPrior #5 

