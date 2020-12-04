% Correction of psi angle based on filament polarity estimation

%% Mean psi angle as prior and removal of flipped segments

% Construct base path
actin_polarity_construct_base_path;

% Define extracted filaments
for t=1:7
    FilamentStruct{t} = ['/home/Medalia/Projects7/Bruno/ActinPolarity/mapping/mapping_144/filaments/filaments_144_tomo_' num2str(t) '.mat'];
end

% Process bundles

zaehler_particles = 1;

group_indx = 1;

for k=[1,2,3,4,5,6,7]
    
    % Load refined filaments
    load(FilamentStruct{k},'filament_struct_ref');
    
    % Process filaments of that bundle
    for i=1:size(filament_struct_ref,2)
    
         % Extract filament spin
         spin_filament_red_red = filament_struct_ref(i).spin_filament_red_red;
         
         % Extract filament spin vector
         spin_filament = filament_struct_ref(i).spin_filament;
         
         % Extract segments with psi angles pointing in the most likely direction
         if spin_filament_red_red == 1
              indx_segments = find(spin_filament == 1);
         elseif spin_filament_red_red == 0
              indx_segments = find(spin_filament == 0);
         end
         
         % Calculate the mean psi of those psi angles
         mean_psi_prior = mean(filament_struct_ref(i).actin_data_matrix(indx_segments,29),1);
                  
         % Process those segments
         for j=1:size(indx_segments,1)
            
              % Extract particle_t
              particle_t = filament_struct_ref(i).actin_data_matrix(indx_segments(j),1);
              
              % Extract particle_p
              particle_p = filament_struct_ref(i).actin_data_matrix(indx_segments(j),2);
              
              % Extract particle_n
              particle_n = filament_struct_ref(i).actin_data_matrix(indx_segments(j),3);
              
              % Create particle name
              particle_name = actin_polarity_construct_particle_identifier(particle_t,particle_p,particle_n);
              
              % Extract particle index
              particle_name_pre = particle_name(regexp(particle_name,'_n_')+3:end);
              if size(particle_name_pre,2) == 5
                  particle_name_rln = ['0' particle_name_pre];
              elseif size(particle_name_pre,2) == 6
                  particle_name_rln = particle_name_pre;
              end
              
              % Extract psi angle
              psi_prior_rev = filament_struct_ref(i).actin_data_matrix(indx_segments(j),29);
              
              % Prepare Relion starfile
              
              if group_indx == 1
              
                  rln_starfile{zaehler_particles} = [particle_name_rln '@/home/Medalia/Projects7/Bruno/ActinScratchX03/ActinPolarity/particles/actin_polarity_particles_list_bin0_144_proj_norm_inv.mrcs' '   ' ...
                                                                    'actin_polarity_tomo_' num2str(group_indx) '   ' ...
                                                 num2str(round((rand(1)-0.5).*360)) '   ' ...
                                                 num2str(90) '   ' ...
                                                 num2str(mean_psi_prior)];%num2str(psi_prior_rev)
              group_indx = 2;
              
              else
                  
                  rln_starfile{zaehler_particles} = [particle_name_rln '@/home/Medalia/Projects7/Bruno/ActinScratchX03/ActinPolarity/particles/actin_polarity_particles_list_bin0_144_proj_norm_inv.mrcs' '   ' ...
                                                                    'actin_polarity_tomo_' num2str(group_indx) '   ' ...
                                                 num2str(round((rand(1)-0.5).*360)) '   ' ...
                                                 num2str(90) '   ' ...
                                                 num2str(mean_psi_prior)];%num2str(psi_prior_rev)
              group_indx = 1;
              
              end
              
              zaehler_particles = zaehler_particles + 1;
              
         end
        
    end
    
    disp(k);
    
end

% % Write starfile
% placeh = num2str(0);
%        save('/home/Medalia/Projects7/Bruno/ActinPolarity/mapping/mapping_144/psi_prior/particles_job004_psi_prior_rev.star', 'placeh', '-ASCII');
% fid = fopen('/home/Medalia/Projects7/Bruno/ActinPolarity/mapping/mapping_144/psi_prior/particles_job004_psi_prior_rev.star', 'w');
% for k=1:size(rln_starfile,2)
%     fprintf(fid, '%s\n', char(rln_starfile{k}));
% end
% fclose(fid);

% Add header manually

% data_
%
% loop_
% _rlnImageName #1
% _rlnMicrographName #2
% _rlnAngleRotPrior #3
% _rlnAngleTiltPrior #4
% _rlnAnglePsiPrior #5


%% Mean psi angle as prior and flipping of segments

% Construct base path
actin_polarity_construct_base_path;

% Define extracted filaments
for t=1:7
    FilamentStruct{t} = ['/home/Medalia/Projects7/Bruno/ActinPolarity/mapping/mapping_144/filaments/filaments_144_tomo_' num2str(t) '.mat'];
end

% Process bundles

zaehler_particles = 1;

group_indx = 1;

for k=[1,2,3,4,5,6,7]
    
    
    % Load refined filaments
    load(FilamentStruct{k},'filament_struct_ref');
    
    % Process filaments of that bundle
    for i=1:size(filament_struct_ref,2)
    
         % Extract filament spin
         spin_filament_red_red = filament_struct_ref(i).spin_filament_red_red;
         
         % Extract filament spin vector
         spin_filament = filament_struct_ref(i).spin_filament;
         
         % Extract segments with psi angles pointing in the most likely direction
         if spin_filament_red_red == 1
              indx_segments = find(spin_filament == 1);
         elseif spin_filament_red_red == 0
              indx_segments = find(spin_filament == 0);
         end
         
         % Calculate the mean psi of those psi angles
         mean_psi_prior = mean(filament_struct_ref(i).actin_data_matrix(indx_segments,29),1);
                  
         % Process those segments
         for j=1:size(spin_filament,1)
            
              % Extract particle_t
              particle_t = filament_struct_ref(i).actin_data_matrix(j,1);
              
              % Extract particle_p
              particle_p = filament_struct_ref(i).actin_data_matrix(j,2);
              
              % Extract particle_n
              particle_n = filament_struct_ref(i).actin_data_matrix(j,3);
              
              % Create particle name
              particle_name = actin_polarity_construct_particle_identifier(particle_t,particle_p,particle_n);
              
              % Extract particle index
              particle_name_pre = particle_name(regexp(particle_name,'_n_')+3:end);
              if size(particle_name_pre,2) == 5
                  particle_name_rln = ['0' particle_name_pre];
              elseif size(particle_name_pre,2) == 6
                  particle_name_rln = particle_name_pre;
              end
              
              % Extract psi angle
              psi_prior_rev = filament_struct_ref(i).actin_data_matrix(j,29);
              
              % Prepare Relion starfile
              
              if group_indx == 1
              
                  rln_starfile{zaehler_particles} = [particle_name_rln '@/home/Medalia/Projects7/Bruno/ActinScratchX03/ActinPolarity/particles/actin_polarity_particles_list_bin0_144_proj_norm_inv.mrcs' '   ' ...
                                                                    'actin_polarity_tomo_' num2str(group_indx) '   ' ...
                                                 num2str(round((rand(1)-0.5).*360)) '   ' ...
                                                 num2str(90) '   ' ...
                                                 num2str(mean_psi_prior)];%num2str(ps/mapping_144/psi_prior/particles_job004_psi_prior_flii_prior_rev)
              group_indx = 2;
              
              else
                  
                  rln_starfile{zaehler_particles} = [particle_name_rln '@/home/Medalia/Projects7/Bruno/ActinScratchX03/ActinPolarity/particles/actin_polarity_particles_list_bin0_144_proj_norm_inv.mrcs' '   ' ...
                                                                    'actin_polarity_tomo_' num2str(group_indx) '   ' ...
                                                 num2str(round((rand(1)-0.5).*360)) '   ' ...
                                                 num2str(90) '   ' ...
                                                 num2str(mean_psi_prior)];%num2str(psi_prior_rev)
              group_indx = 1;
              
              end
              
              zaehler_particles = zaehler_particles + 1;
              
         end
        
    end
    
    disp(k);
    
end

% Write starfile
placeh = num2str(0);
       save('/home/Medalia/Projects7/Bruno/ActinPolarity/mapping/mapping_144/psi_prior/particles_job004_psi_prior_flip_rev.star', 'placeh', '-ASCII');
fid = fopen('/home/Medalia/Projects7/Bruno/ActinPolarity/mapping/mapping_144/psi_prior/particles_job004_psi_prior_flip_rev.star', 'w');
for k=1:size(rln_starfile,2)
    fprintf(fid, '%s\n', char(rln_starfile{k}));
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

