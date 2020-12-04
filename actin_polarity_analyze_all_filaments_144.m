
% Construct base path
actin_polarity_construct_base_path;

% Define extracted filaments
for t=1:7
    FilamentStruct{t} = ['/home/Medalia/Projects7/Bruno/ActinPolarity/mapping/mapping_144/filaments/filaments_144_tomo_' num2str(t) '.mat'];
end

% Combine all bundles
filament_struct_all = [];
for k=[1,2,3,4,5,6,7]

    % Load filament structure
    load(FilamentStruct{k},'filament_struct_ref');
    
    % Add tomogram index
    for i=1:size(filament_struct_ref,2)
         filament_struct_ref(i).tomogram_index = k;
    end
    
    filament_struct_all = [filament_struct_all filament_struct_ref];
    
end

% Extract length of selected filaments

all_filament_length_sel = [];

for k=1:size(filament_struct_all,2)

    all_filament_length_sel = [all_filament_length_sel filament_struct_all(k).length_filament];
    
end

% Convert to nm units
all_filament_length_sel_nm = (4.*3.443./10).*all_filament_length_sel;

% Plot histogram
figure;figure_indx = gcf;figure(figure_indx);histogram(all_filament_length_sel_nm);
xlabel('Length of rec. filaments (nm)');
ylabel('Number of rec. filaments');

% Extract distance between segments of selected filaments

all_filament_seg_dist_sel = [];

for k=1:size(filament_struct_all,2)

    all_filament_seg_dist_sel = [all_filament_seg_dist_sel filament_struct_all(k).length_filament_vect'];
    
end

% Convert to nm units
all_filament_seg_dist_sel_nm = (4.*3.443./10).*all_filament_seg_dist_sel;

% Plot histogram
figure;figure_indx = gcf;figure(figure_indx);histogram(all_filament_seg_dist_sel_nm);
xlim([5 60]);
xticks([11 22 33 44 55]);
xlabel('Distance between segments (nm)');
ylabel('Number of segments');

% Extract psi angles and spins of selected filaments

xxx = 0;
yyy = 0;

filament_struct_psi_spin_sel = [];
for k=[1,2,3,4,5,6,7]

    % Load filament structure
    load(FilamentStruct{k},'filament_struct_ref');
    
    % Initialize selected psi angles
    psi_angles_sel = [];
    
    % Initialize reconstructed and selected filament spins
    filament_spin_sel = [];
    
    % Extract psi angles of selected filaments
    for i=1:size(filament_struct_ref,2)
         
         psi_angles_sel = [psi_angles_sel filament_struct_ref(i).psi_filament'];
              
         filament_spin_sel = [filament_spin_sel filament_struct_ref(i).spin_filament_red_red];
         
    end
    
    filament_struct_psi_spin_sel(k).psi_angles_sel = psi_angles_sel;
    clear psi_angles_sel;
    
    filament_struct_psi_spin_sel(k).filament_spin_sel = filament_spin_sel;
    clear filament_spin_sel;
    
    disp('----------------');
    disp(['Tomogram index: ' num2str(k)]);
    disp(['Number of rec. filaments with spin 0: ' num2str(size(find(filament_struct_psi_spin_sel(k).filament_spin_sel==0),2))]);
    disp(['Percentage of rec. filaments with spin 0: ' num2str(round(100.*(size(find(filament_struct_psi_spin_sel(k).filament_spin_sel==0),2)./(size(find(filament_struct_psi_spin_sel(k).filament_spin_sel==0),2)+size(find(filament_struct_psi_spin_sel(k).filament_spin_sel==1),2)))))]);
    disp(['Number of rec. filaments with spin 1: ' num2str(size(find(filament_struct_psi_spin_sel(k).filament_spin_sel==1),2))]);
    disp(['Percentage of rec. filaments with spin 1: ' num2str(round(100.*(size(find(filament_struct_psi_spin_sel(k).filament_spin_sel==1),2)./(size(find(filament_struct_psi_spin_sel(k).filament_spin_sel==0),2)+size(find(filament_struct_psi_spin_sel(k).filament_spin_sel==1),2)))))]);
    disp('----------------');
    
    xxx = xxx + size(filament_struct_psi_spin_sel(k).psi_angles_sel,2);
    yyy = yyy + size(filament_struct_psi_spin_sel(k).filament_spin_sel,2);
    
end

