
% Define extracted filaments
for t=1:7
    FilamentStruct{t} = ['/home/Medalia/Projects7/Bruno/ActinPolarity/mapping/mapping_144/filaments/filaments_144_tomo_' num2str(t) '.mat'];
end

% Initialize all distances of filaments from the central segment
all_neighbourhood_distances = [];

% Initialize all number of filaments within the neighbourhood
all_filaments_within_neighbourhood = [];

% Process bundles
spin_neighbourhood_matrix = [];
zaehler = 1;
for k=[1,2,3,4,5,6,7]

    % Load filament structure
    load(FilamentStruct{k},'filament_struct_ref');
    
    % Label spin neighbourhoods
    for i=1:size(filament_struct_ref,2)
      for j=1:size(filament_struct_ref(i).segment_neighbourhood_epsilon,2)
          
          % Extract spin neighbourhood
          spin_neighbourhood = filament_struct_ref(i).segment_neighbourhood_epsilon(j).spin_neighbourhood;
          
          % Extract distance neighbourhood
          dist_neighbourhood = filament_struct_ref(i).segment_neighbourhood_epsilon(j).dist_neighbourhood;
          
          % Extract the spin of the central segment of the neighbourhood
          segment_spin = filament_struct_ref(i).spin_filament_red_vect(j,1);
          
          if isempty(spin_neighbourhood)
              
              spin_neighbourhood_matrix(zaehler,1) = -1;% Mean spin of filament neighbourhood
              spin_neighbourhood_matrix(zaehler,2) = -1;% Number of filaments within probed neighbourhood
              spin_neighbourhood_matrix(zaehler,3) = -1;% Mean distance of filaments from the central segment in the probed neighbourhood
              spin_neighbourhood_matrix(zaehler,4) = -1;% Std of distances of the filaments from the central segment in the probed neighbourhood
              spin_neighbourhood_matrix(zaehler,5) = -1;% Central segment spin
              spin_neighbourhood_matrix(zaehler,6) = -1;% Neighbourhood spin compared to central segment spin
              
              spin_neighbourhood_matrix(zaehler,7) = k;% Bundle identifier
              spin_neighbourhood_matrix(zaehler,8) = i;% Filament identifier
              spin_neighbourhood_matrix(zaehler,9) = j;% Segment identifier
              
              zaehler = zaehler + 1;
              
          else
              
              spin_neighbourhood_matrix(zaehler,1) = mean(spin_neighbourhood);                % Mean spin of filament neighbourhood
              spin_neighbourhood_matrix(zaehler,2) = size(spin_neighbourhood,2);              % Number of filaments within probed neighbourhood
              spin_neighbourhood_matrix(zaehler,3) = mean(dist_neighbourhood);                % Mean distance of filaments from the central segment in the probed neighbourhood
              spin_neighbourhood_matrix(zaehler,4) = std(dist_neighbourhood);                 % Std of distances of the filaments from the central segment in the probed neighbourhood
              spin_neighbourhood_matrix(zaehler,5) = segment_spin;                            % Central segment spin
              spin_neighbourhood_matrix(zaehler,6) = mean(spin_neighbourhood) - segment_spin; % Neighbourhood spin compared to central segment spin
              
              spin_neighbourhood_matrix(zaehler,7) = k;% Bundle identifier
              spin_neighbourhood_matrix(zaehler,8) = i;% Filament identifier
              spin_neighbourhood_matrix(zaehler,9) = j;% Segment identifier
              
              % Collect all distances of filaments from the central segment
              all_neighbourhood_distances = [all_neighbourhood_distances dist_neighbourhood];
              
              % Collect all number of filaments within the neighbourhood
              all_filaments_within_neighbourhood = [all_filaments_within_neighbourhood size(spin_neighbourhood,2).*ones(1,size(dist_neighbourhood,2))];
              
              zaehler = zaehler + 1;
              
          end
          
      end
    end
    
    disp(k);
    
end

% % Plot all neighbourhoods distance distribution
% figure;figure_indx = gcf;figure(figure_indx);histogram(all_neighbourhood_distances);
% xlabel('Filament distance in neighbourhood (pixels)');
% xlim([0 30]);
% ylim([0 2200]);

% % Plot all number of filaments per neighbourhood distribution
% figure;figure_indx = gcf;figure(figure_indx);histogram(all_filaments_within_neighbourhood);
% xlabel('Number of filaments in neigbourhood');
% xlim([0 25]);
% ylim([0 7000]);

% Define maximum number of filaments within neighbourhood
max_number_of_filaments = max(unique(all_filaments_within_neighbourhood));

% Initialize neighbourhood evaluation matrix
neighbourhood_matrix_eval = zeros(7,4);

% Initialize all neighbourhood mean and std of distances of filaments from the segment
all_neighbourhood_mean_distances = [];
all_neighbourhood_std_distances = [];

% Define histogram bins
xxx = -1:0.1:+1; %spin bins
yyy = 1:1:max_number_of_filaments; % number of filaments within neighbourhood

% Perform neighbourhood analysis for each bundle
zaehler = 1;
for k=[1,2,3,4,5,6,7]

    % Select segments of bundle
    indx_segments_bundle = find(spin_neighbourhood_matrix(:,7) == k);
    
    % Reduce neighbourhood matrix to bundle
    spin_neighbourhood_matrix_red = spin_neighbourhood_matrix(indx_segments_bundle,:);
    
    % Select neighbourhoods 
    indx_segments = find(spin_neighbourhood_matrix_red(:,1) >= 0 & spin_neighbourhood_matrix_red(:,2) <= max_number_of_filaments);
    
    % Collect all neighbourhood mean and std of distances
    all_neighbourhood_mean_distances = [all_neighbourhood_mean_distances; spin_neighbourhood_matrix_red(indx_segments,3)];
    all_neighbourhood_std_distances = [all_neighbourhood_std_distances; spin_neighbourhood_matrix_red(indx_segments,4)];
    
    % Calculate mean segment to filament distances per bundle
    mean_segment_to_filament_distance = mean(spin_neighbourhood_matrix_red(indx_segments,3));
    
    % Map spin neighbourhood data to histogram
    bundle_neighbourhood_hist = hist3([spin_neighbourhood_matrix_red(indx_segments,6) spin_neighbourhood_matrix_red(indx_segments,2)],{xxx yyy});
    
    % Plot bundle histograms
    figure(1);
    subplot(1,7,k);imagesc(bundle_neighbourhood_hist);
    axis off;
    axis on;
    set(gca,'YTick',[1 6 11 16 21]);
    set(gca,'YTickLabel',{'-1'         '-0.5'            '0'          '+0.5'            '+1'});
    title(k);
    
    % Sum spin values of different neighbourhoods
    A = sum(bundle_neighbourhood_hist,2)';
    
    % Plot bundle bar histograms
    figure(2);
    subplot(1,7,k);
    bar(A);
    xlim([0 22]);
    set(gca,'XTick',[1 6 11 16 21]);
    set(gca,'XTickLabel',{'-1'         '-0.5'            '0'          '+0.5'            '+1'});
    title(k);
    
    % Calculate number of segments in total parallel or antiparallel neighbourhoods
    neighbourhood_matrix_eval(k,1) = A(1,11);
    neighbourhood_matrix_eval(k,2) = A(1,1) + A(1,21);
    
    % Calculate weighted number of segments per neighbourhood type and bundle
    neighbourhood_matrix_eval(k,3) = round(sum(abs(abs(xxx)-1).*A));
    neighbourhood_matrix_eval(k,4) = round(sum(abs(xxx).*A));
    
    % Calculate bundle neighbourhood scores
    bundle_neighbourhood_scores(1,zaehler) = k;
    
    bundle_neighbourhood_scores(2,zaehler) = neighbourhood_matrix_eval(k,1);    %Number of segments in total parallel neighbourhoods
    bundle_neighbourhood_scores(3,zaehler) = neighbourhood_matrix_eval(k,2);    %Number of segments in total anti-parallel neighbourhoods
                                                                                
    bundle_neighbourhood_scores(4,zaehler) = neighbourhood_matrix_eval(k,3);    %Bundle score parallel
    bundle_neighbourhood_scores(5,zaehler) = neighbourhood_matrix_eval(k,4);    %Bundle score anti-parallel
    
    bundle_neighbourhood_scores(6,zaehler) = round((   neighbourhood_matrix_eval(k,3) - neighbourhood_matrix_eval(k,4)   )./ ...
                                                   (   neighbourhood_matrix_eval(k,3) + neighbourhood_matrix_eval(k,4)   ).*1000)./1000;
    
    
    % Print neighbourhood analysis protocol
    disp('****************************************************************');
    disp(['Bundle index ' num2str(k)]);
    disp('****************************************************************');
    disp(['Number of segments in that bundle: ' num2str(size(spin_neighbourhood_matrix_red,1))]);
    disp(['Number of segments in that bundle with neighbourhood: ' num2str(size(find(spin_neighbourhood_matrix_red(:,1)>=0),1))]);
    disp(['Number of segments in that bundle without neighbourhood: ' num2str(size(spin_neighbourhood_matrix_red,1)-size(find(spin_neighbourhood_matrix_red(:,1)>=0),1))]);
    disp(['Fraction of segments with neighbourhood: ' num2str(size(find(spin_neighbourhood_matrix_red(:,1)>=0),1)./size(spin_neighbourhood_matrix_red,1))]);
    disp(['Number of segments analyzed: ' num2str(size(spin_neighbourhood_matrix_red(indx_segments,:),1))]);
    disp(['Fraction of segments analyzed: ' num2str(size(spin_neighbourhood_matrix_red(indx_segments,:),1)./size(spin_neighbourhood_matrix_red,1))]);
    disp(['Mean segment to filament distance: ' num2str(mean_segment_to_filament_distance) ' pixel / ' num2str(mean_segment_to_filament_distance.*4.*0.3443) ' nm']);
    disp(['Number of segments in total parallel neighbourhoods: ' num2str(neighbourhood_matrix_eval(k,1))]);
    disp(['Number of segments in total anti-parallel neighbourhoods: ' num2str(neighbourhood_matrix_eval(k,2))]);
    disp(['Bundle score parallel: ' num2str(neighbourhood_matrix_eval(k,3))]);
    disp(['Bundle score anti-parallel: ' num2str(neighbourhood_matrix_eval(k,4))]);
    disp(['Bundle topology score: ' num2str(bundle_neighbourhood_scores(6,zaehler))]);
    
    zaehler = zaehler + 1;
    
end

% % Plot all filament mean distances
% figure;figure_indx = gcf;figure(figure_indx);histogram(all_neighbourhood_mean_distances.*4.*0.3443);
% xlabel('Mean filament distance in neighbourhood (nm)');
% xlim([0 40]);
% ylim([0 1100]);

% Sort according topology score
[~,indx_sort_topology_score] = sort(bundle_neighbourhood_scores(6,:));
bundle_neighbourhood_scores_sort = bundle_neighbourhood_scores(:,indx_sort_topology_score);

% Display sorted topology score and sorted tomogram indices
bundle_neighbourhood_scores_sort(6,:)
bundle_neighbourhood_scores_sort(1,:)

