%% Create barbed end Volume

% Load actin average
actin_avg = tom_mrcread('/home/Medalia/Projects7/Bruno/ActinPolarity/averaging/rln_144/Refine3D/job006/run_class001.mrc');
actin_avg = actin_avg.Value;

% Create the mask
mask = ones(144,144,144);
mask = single(mask);

% Define the mask
mask(:,:,25:144) = 0;

% Multiply mask with actin_avg
merge = actin_avg.*mask;

% Save masked volume
tom_mrcwrite(merge,'name','/home/Medalia/Projects7/Bruno/ActinPolarity/mapping/mapping_144/plus_end_assignment/rsc/merge.mrc');


%% Load actin data matrix

% Define psi thresholds
psi_threshold = zeros(7,6);

psi_threshold(1,1)  = -73;   psi_threshold(1,2)  = -17;   psi_threshold(1,3)  = +122;   psi_threshold(1,4)  = +158;    psi_threshold(1,5) = -73;   psi_threshold(1,6)  = -17;

psi_threshold(2,1)  = -180;  psi_threshold(2,2)  = -149;  psi_threshold(2,3)  = -20;    psi_threshold(2,4)  = +28;    psi_threshold(2,5)  = +146;  psi_threshold(2,6)  = +180;

psi_threshold(3,1)  = -141;  psi_threshold(3,2)  = -107;  psi_threshold(3,3)  = +40;    psi_threshold(3,4)  = +77;    psi_threshold(3,5)  = -141;  psi_threshold(3,6)  = -107;

psi_threshold(4,1)  = -148;  psi_threshold(4,2)  = -96;   psi_threshold(4,3)  = +24;    psi_threshold(4,4)  = +85;    psi_threshold(4,5)  = -148;  psi_threshold(4,6)  = -96;

psi_threshold(5,1)  = -132;  psi_threshold(5,2)  = -97;   psi_threshold(5,3)  = +37;    psi_threshold(5,4)  = +84;    psi_threshold(5,5)  = -132;  psi_threshold(5,6)  = -97;

psi_threshold(6,1)  = -83;   psi_threshold(6,2)  = -36;   psi_threshold(6,3)  = +101;   psi_threshold(6,4)  = +143;   psi_threshold(6,5)  = -83;   psi_threshold(6,6)  = -36;

psi_threshold(7,1)  = -76;   psi_threshold(7,2)  = -26;   psi_threshold(7,3)  = +103;   psi_threshold(7,4)  = +155;   psi_threshold(7,5)  = -76;   psi_threshold(7,6)  = -26;

% Construct base path
actin_polarity_construct_base_path;

% Load actin data matrix (created in "actin_polarity_generate_data_matrix_144.m")
load('/home/Medalia/Projects7/Bruno/ActinPolarity/mapping/mapping_144/actin_polarity_mapping_144_step_2.mat','actin_data_matrix');

% Extract bundles
indx_bundle_n = unique(actin_data_matrix(:,1))'; % unique indices tomogram_n

% Plot actin plot
for k=2:2

disp (k);  

% Find indices of bundle that is analyzed minus excluded filaments
indx_bundle_x = find(actin_data_matrix(:,1)' == indx_bundle_n(1,k) & actin_data_matrix(:,4)' == 1);

% Extract angles
psi_angle = actin_data_matrix(indx_bundle_x,19);

% Create chimera volumes - job006

% Load actin average
actin_avg = tom_mrcread('/home/Medalia/Projects7/Bruno/ActinPolarity/averaging/rln_144/Refine3D/job006/run_class001.mrc');
actin_avg = tom_bin(actin_avg.Value,3);

% Initialize volume / red filaments
vol_a = zeros(512,512,512);
vol_b = zeros(512,512,512);

for i=1:size(indx_bundle_x,2)
    
    % Inversly align average
    particle_align = double(tom_rotate(tom_shift(actin_avg,[actin_data_matrix(indx_bundle_x(i),21)./8 actin_data_matrix(indx_bundle_x(i),22)./8 actin_data_matrix(indx_bundle_x(i),23)./8]),[actin_data_matrix(indx_bundle_x(i),18) actin_data_matrix(indx_bundle_x(i),19) actin_data_matrix(indx_bundle_x(i),20)]));

    % Normalize aligned particle
    particle_align = (particle_align - mean(particle_align(:)))./std((particle_align(:)));

    % Initialize aligned particle volume
    vol_particle_aligned = zeros(512,512,512);    
    
     if     psi_angle(i) >= psi_threshold(k,1) && psi_angle(i) <= psi_threshold(k,2)
         try
              vol_particle_aligned(round(actin_data_matrix(indx_bundle_x(i),8)./8)-9+1:round(actin_data_matrix(indx_bundle_x(i),8)./8)+9,round(actin_data_matrix(indx_bundle_x(i),9)./8)-9+1:round(actin_data_matrix(indx_bundle_x(i),9)./8)+9,round(actin_data_matrix(indx_bundle_x(i),10)./8)-9+1:round(actin_data_matrix(indx_bundle_x(i),10)./8)+9) = particle_align;
         end
         vol_a = vol_a + vol_particle_aligned;%spin=0/blue
                    
    elseif psi_angle(i) >= psi_threshold(k,3) && psi_angle(i) <= psi_threshold(k,4)
         try
              vol_particle_aligned(round(actin_data_matrix(indx_bundle_x(i),8)./8)-9+1:round(actin_data_matrix(indx_bundle_x(i),8)./8)+9,round(actin_data_matrix(indx_bundle_x(i),9)./8)-9+1:round(actin_data_matrix(indx_bundle_x(i),9)./8)+9,round(actin_data_matrix(indx_bundle_x(i),10)./8)-9+1:round(actin_data_matrix(indx_bundle_x(i),10)./8)+9) = particle_align;
         end
         vol_b = vol_b + vol_particle_aligned;%spin=1/red
                    
    elseif psi_angle(i) >= psi_threshold(k,5) && psi_angle(i) <= psi_threshold(k,6)
                    
         try
              vol_particle_aligned(round(actin_data_matrix(indx_bundle_x(i),8)./8)-9+1:round(actin_data_matrix(indx_bundle_x(i),8)./8)+9,round(actin_data_matrix(indx_bundle_x(i),9)./8)-9+1:round(actin_data_matrix(indx_bundle_x(i),9)./8)+9,round(actin_data_matrix(indx_bundle_x(i),10)./8)-9+1:round(actin_data_matrix(indx_bundle_x(i),10)./8)+9) = particle_align;
         end
         vol_a = vol_a + vol_particle_aligned;%spin=0/blue

    end
             
    disp(i);
end

% Save volume
tom_emwritec3(['/home/Medalia/Projects7/Bruno/ActinPolarity/mapping/mapping_144/plus_end_assignment/actin_avg_a_' num2str(k) '.em'],single(vol_a));
tom_emwritec3(['/home/Medalia/Projects7/Bruno/ActinPolarity/mapping/mapping_144/plus_end_assignment/actin_avg_b_' num2str(k) '.em'],single(vol_b));


end


%% Load actin data matrix

% Define psi thresholds
psi_threshold = zeros(7,6);

psi_threshold(1,1)  = -73;   psi_threshold(1,2)  = -17;   psi_threshold(1,3)  = +122;   psi_threshold(1,4)  = +158;    psi_threshold(1,5) = -73;   psi_threshold(1,6)  = -17;

psi_threshold(2,1)  = -180;  psi_threshold(2,2)  = -149;  psi_threshold(2,3)  = -20;    psi_threshold(2,4)  = +28;    psi_threshold(2,5)  = +146;  psi_threshold(2,6)  = +180;

psi_threshold(3,1)  = -141;  psi_threshold(3,2)  = -107;  psi_threshold(3,3)  = +40;    psi_threshold(3,4)  = +77;    psi_threshold(3,5)  = -141;  psi_threshold(3,6)  = -107;

psi_threshold(4,1)  = -148;  psi_threshold(4,2)  = -96;   psi_threshold(4,3)  = +24;    psi_threshold(4,4)  = +85;    psi_threshold(4,5)  = -148;  psi_threshold(4,6)  = -96;

psi_threshold(5,1)  = -132;  psi_threshold(5,2)  = -97;   psi_threshold(5,3)  = +37;    psi_threshold(5,4)  = +84;    psi_threshold(5,5)  = -132;  psi_threshold(5,6)  = -97;

psi_threshold(6,1)  = -83;   psi_threshold(6,2)  = -36;   psi_threshold(6,3)  = +101;   psi_threshold(6,4)  = +143;   psi_threshold(6,5)  = -83;   psi_threshold(6,6)  = -36;

psi_threshold(7,1)  = -76;   psi_threshold(7,2)  = -26;   psi_threshold(7,3)  = +103;   psi_threshold(7,4)  = +155;   psi_threshold(7,5)  = -76;   psi_threshold(7,6)  = -26;

% Construct base path
actin_polarity_construct_base_path;

% Load actin data matrix (created in "actin_polarity_generate_data_matrix_144.m")
load('/home/Medalia/Projects7/Bruno/ActinPolarity/mapping/mapping_144/actin_polarity_mapping_144_step_2.mat','actin_data_matrix');

% Extract bundles
indx_bundle_n = unique(actin_data_matrix(:,1))'; % unique indices tomogram_n

% Plot actin plot
for k=2:2

disp (k);  

% Find indices of bundle that is analyzed minus excluded filaments
indx_bundle_x = find(actin_data_matrix(:,1)' == indx_bundle_n(1,k) & actin_data_matrix(:,4)' == 1);

% Extract angles
psi_angle = actin_data_matrix(indx_bundle_x,19);

% Create chimera volumes - job006

% Load actin average
actin_avg = tom_mrcread('/home/Medalia/Projects7/Bruno/ActinPolarity/mapping/mapping_144/plus_end_assignment/rsc/merge_barbed.mrc');
actin_avg = tom_bin(actin_avg.Value,3);

% Initialize volume / red filaments
vol_a = zeros(512,512,512);
vol_b = zeros(512,512,512);

for i=1:size(indx_bundle_x,2)
    
    % Inversly align average
    particle_align = double(tom_rotate(tom_shift(actin_avg,[actin_data_matrix(indx_bundle_x(i),21)./8 actin_data_matrix(indx_bundle_x(i),22)./8 actin_data_matrix(indx_bundle_x(i),23)./8]),[actin_data_matrix(indx_bundle_x(i),18) actin_data_matrix(indx_bundle_x(i),19) actin_data_matrix(indx_bundle_x(i),20)]));

    % Normalize aligned particle
    particle_align = (particle_align - mean(particle_align(:)))./std((particle_align(:)));

    % Initialize aligned particle volume
    vol_particle_aligned = zeros(512,512,512);    
    
     if     psi_angle(i) >= psi_threshold(k,1) && psi_angle(i) <= psi_threshold(k,2)
         try
              vol_particle_aligned(round(actin_data_matrix(indx_bundle_x(i),8)./8)-9+1:round(actin_data_matrix(indx_bundle_x(i),8)./8)+9,round(actin_data_matrix(indx_bundle_x(i),9)./8)-9+1:round(actin_data_matrix(indx_bundle_x(i),9)./8)+9,round(actin_data_matrix(indx_bundle_x(i),10)./8)-9+1:round(actin_data_matrix(indx_bundle_x(i),10)./8)+9) = particle_align;
         end
         vol_a = vol_a + vol_particle_aligned;
                    
    elseif psi_angle(i) >= psi_threshold(k,3) && psi_angle(i) <= psi_threshold(k,4)
         try
              vol_particle_aligned(round(actin_data_matrix(indx_bundle_x(i),8)./8)-9+1:round(actin_data_matrix(indx_bundle_x(i),8)./8)+9,round(actin_data_matrix(indx_bundle_x(i),9)./8)-9+1:round(actin_data_matrix(indx_bundle_x(i),9)./8)+9,round(actin_data_matrix(indx_bundle_x(i),10)./8)-9+1:round(actin_data_matrix(indx_bundle_x(i),10)./8)+9) = particle_align;
         end
         vol_b = vol_b + vol_particle_aligned;
                    
    elseif psi_angle(i) >= psi_threshold(k,5) && psi_angle(i) <= psi_threshold(k,6)
                    
         try
              vol_particle_aligned(round(actin_data_matrix(indx_bundle_x(i),8)./8)-9+1:round(actin_data_matrix(indx_bundle_x(i),8)./8)+9,round(actin_data_matrix(indx_bundle_x(i),9)./8)-9+1:round(actin_data_matrix(indx_bundle_x(i),9)./8)+9,round(actin_data_matrix(indx_bundle_x(i),10)./8)-9+1:round(actin_data_matrix(indx_bundle_x(i),10)./8)+9) = particle_align;
         end
         vol_a = vol_a + vol_particle_aligned;

    end
             
    disp(i);
end

% Save volume
tom_emwritec3(['/home/Medalia/Projects7/Bruno/ActinPolarity/mapping/mapping_144/plus_end_assignment/masked_a_' num2str(k) '.em'],single(vol_a));
tom_emwritec3(['/home/Medalia/Projects7/Bruno/ActinPolarity/mapping/mapping_144/plus_end_assignment/masked_b_' num2str(k) '.em'],single(vol_b));


end

