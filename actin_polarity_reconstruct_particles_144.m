
% Generate BasePath
actin_polarity_construct_base_path;

% Define microscope parameter
Objectpixelsize = 0.3443;%nm

% Define particle path / bin0_144
ParticlesPath = '/home/Medalia/scratchX03/ActinPolarity/particles/actin_particles_bin0_144/';

% Perform particle reconstruction
for k=1:size(BasePath,2)
    
    % Check if particles are selected
    try
         load([BasePath{k} '/cor/actin_polarity_cor_final_144.mat']);
         particles_selected = 1;
         clear plist;
    catch
         particles_selected = 0;
    end
    
    % Perform reconstruction of selected particles
    if particles_selected == 1
          actin_polarity_perform_particle_rec_144(BasePath{k},ParticlesPath,['actin_particle_t_' num2str(k) '_p_']);
    end
    
    clear particles_selected;
    
    disp(k);
    
end

