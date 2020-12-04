
% Load actin data matrix
load('/home/Medalia/Projects7/Bruno/ActinPolarity/mapping/mapping_144/actin_polarity_mapping_144_step_2.mat','actin_data_matrix');

% Define filament length in terms of segments to analyze threshold
filament_length_thres = 3;

% Define combined confidence score threshold
filament_com_conf_thres = -1;%0,0.6

% Define maximum likelihood estimate confidence interval threshold
filament_mle_conf_int_a_thres = 0.5;%0

% Define radius of epsilon-sphere
epsilon_sphere_radius = 30;%~40 nm radius

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

% Define extracted filaments
for t=1:7
    FilamentStruct{t} = ['/home/Medalia/Projects7/Bruno/ActinPolarity/mapping/mapping_144/filaments/filaments_144_tomo_' num2str(t) '.mat'];
end

% Initialize overall distributions and tomo membership
psi_distribution_all = [];
ccc_distribution_all = [];
segment_sensitivity_distribution_all = [];
tomo_distribution_all = [];

% Initialize all filament confidence values
spin_filament_conf_all = [];
spin_filament_conf_uni_all = [];
spin_filament_conf_all_0 = [];
spin_filament_conf_all_1 = [];
spin_filament_conf_seg_sense_all = [];
spin_filament_conf_mean_ccc_all = [];
spin_filament_conf_com_score_all = [];
spin_filament_mle_conf_int_a_all = [];
spin_filament_conf_tomo_all = [];

% Number of segments inside/outside of distribution
zaehler_seg_in = 0;
zaehler_seg_out = 0;

% Process bundles
for k=[1,2,3,4,5,6,7]

    % Select only those segments that were used in the respective tomogram and in Relion 3D refine
    indx_segments = find(actin_data_matrix(:,1) == k & actin_data_matrix(:,4) == 1);
    
    % Reduce actin data matrix
    actin_data_matrix_red = actin_data_matrix(indx_segments,:);
    
    % Extract filament assignment of each segment in the respective tomogram
    filament_vect = actin_data_matrix_red(:,38);
    
    % Extract filaments and filament distributions
    filament_struct = [];
    psi_distribution = [];
    ccc_distribution = [];
    tomo_distribution = [];
    segment_sensitivity_distribution = [];
    zaehler_f = 1;
    
    for i=1:max(filament_vect(:))
        
        % Extract filament indices
        indx_filament_ext = find(filament_vect==i);
        
        if ~isempty(indx_filament_ext) && size(indx_filament_ext,1) >= filament_length_thres
        
              filament_struct(zaehler_f).filament_indx = zaehler_f;
              
              % Extract coordinates of alignment corrected segments
              cor_filament_ext = [actin_data_matrix_red(indx_filament_ext,8)./4 actin_data_matrix_red(indx_filament_ext,9)./4 actin_data_matrix_red(indx_filament_ext,10)./4];
              
              % Extract actin directions
              psi_filament_ext = actin_data_matrix_red(indx_filament_ext,19);
              psi_distribution = [psi_distribution; psi_filament_ext];
              
              % Extract ccc distribution
              ccc_filament_ext = actin_data_matrix_red(indx_filament_ext,35);
              ccc_distribution = [ccc_distribution; ccc_filament_ext];
              
              % Extract tomo distribution
              tomo_filament_ext = actin_data_matrix_red(indx_filament_ext,1);
              tomo_distribution = [tomo_distribution; tomo_filament_ext];
              
              % Extract filament sensitivity distribution
              segment_sensitivity_filament_ext = actin_data_matrix_red(indx_filament_ext,37);
              segment_sensitivity_distribution = [segment_sensitivity_distribution; segment_sensitivity_filament_ext];
              
              % Extract actin vector
              vect_filament_ext = [actin_data_matrix_red(indx_filament_ext,24) actin_data_matrix_red(indx_filament_ext,25) actin_data_matrix_red(indx_filament_ext,26)];
              
              % Reduce actin data matrix to this filament
              actin_data_matrix_ext = actin_data_matrix_red(indx_filament_ext,:);
              
              filament_struct(zaehler_f).cor_filament_ext = cor_filament_ext;
              filament_struct(zaehler_f).psi_filament_ext = psi_filament_ext;
              filament_struct(zaehler_f).ccc_filament_ext = ccc_filament_ext;
              filament_struct(zaehler_f).segment_sensitivity_filament_ext = segment_sensitivity_filament_ext;
              filament_struct(zaehler_f).vect_filament_ext = vect_filament_ext;
              filament_struct(zaehler_f).actin_data_matrix_ext = actin_data_matrix_ext;
              
              zaehler_f = zaehler_f + 1;
              
        end
        
    end
    
%     % Plot psi/ccc/segment_sensitivity distributions
%     figure;
%     subplot(1,3,1);histogram(psi_distribution(:),50);xlim([-200 +200]);
%     subplot(1,3,2);histogram(ccc_distribution(:),50);xlim([0 +0.3]);
%     subplot(1,3,3);histogram(segment_sensitivity_distribution,50);xlim([-1.2 +1.2]);
%     title(['Actin polarity tomogram: ' num2str(k)]);
    
%     % Plot psi/ccc/segment_sensitivity distributions
%     figure(1);
%     subplot(1,7,k);histogram(psi_distribution(:),50);xlim([-200 +200]);title(['Psi (' num2str(k) ')']);
%     figure(2);
%     subplot(1,7,k);histogram(ccc_distribution(:),50);xlim([0 +0.3]);title(['CCC (' num2str(k) ')']);
%     figure(3);
%     subplot(1,7,k);histogram(segment_sensitivity_distribution,50);xlim([-1.2 +1.2]);title(['Seg-sens (' num2str(k) ')']);

%     % Plot psi distributions
%     figure(k);histogram(psi_distribution(:),90);xlim([-200 +200]);title(['Psi (' num2str(k) ')']);
    
    %continue;
    
    % Collect all distribution values and tomo membership
    psi_distribution_all = [psi_distribution_all; psi_distribution];
    ccc_distribution_all = [ccc_distribution_all; ccc_distribution];
    segment_sensitivity_distribution_all = [segment_sensitivity_distribution_all; segment_sensitivity_distribution];
    tomo_distribution_all = [tomo_distribution_all; tomo_distribution];
    
    % Refine filaments and create spin vectors
    filament_struct_ref = [];
    zaehler_f_ref = 1;
    
    zaehler_seg_tomo_spin_0 = 0;
    zaehler_seg_tomo_spin_1 = 0;
    
    for i=1:size(filament_struct,2)
       
        zaehler_s = 1;
        
        filament_spin_indx_vect = [];
        filament_spin_vect = [];
        
        for j=1:size(filament_struct(i).cor_filament_ext,1)
           
            if     filament_struct(i).psi_filament_ext(j,1) >= psi_threshold(k,1) && filament_struct(i).psi_filament_ext(j,1) <= psi_threshold(k,2)
                    
                        filament_spin_indx_vect(1,zaehler_s) = j;
                        zaehler_s = zaehler_s + 1;
                        filament_spin_vect(1,j) = 0;
                        
                        zaehler_seg_tomo_spin_0 = zaehler_seg_tomo_spin_0 + 1;
                        
                        zaehler_seg_in = zaehler_seg_in + 1;
                    
            elseif filament_struct(i).psi_filament_ext(j,1) >= psi_threshold(k,3) && filament_struct(i).psi_filament_ext(j,1) <= psi_threshold(k,4)
                    
                        filament_spin_indx_vect(1,zaehler_s) = j;
                        zaehler_s = zaehler_s + 1;
                        filament_spin_vect(1,j) = 1;
                        
                        zaehler_seg_tomo_spin_1 = zaehler_seg_tomo_spin_1 + 1;
                        
                        zaehler_seg_in = zaehler_seg_in + 1;
                    
            elseif filament_struct(i).psi_filament_ext(j,1) >= psi_threshold(k,5) && filament_struct(i).psi_filament_ext(j,1) <= psi_threshold(k,6)
                    
                        filament_spin_indx_vect(1,zaehler_s) = j;
                        zaehler_s = zaehler_s + 1;
                        filament_spin_vect(1,j) = 0;
                        
                        zaehler_seg_tomo_spin_0 = zaehler_seg_tomo_spin_0 + 1;
                        
                        zaehler_seg_in = zaehler_seg_in + 1;
                    
            else
                        filament_spin_vect(1,j) = 0.5;%This is the case for ~6% of the segments!
                        
                        zaehler_seg_out = zaehler_seg_out + 1;
                        
            end
            
        end
        
        if size(filament_spin_indx_vect,2) >= filament_length_thres
        
              filament_struct_ref(zaehler_f_ref).cor_filament  = filament_struct(i).cor_filament_ext(filament_spin_indx_vect,:);
              
              filament_struct_ref(zaehler_f_ref).psi_filament  = filament_struct(i).psi_filament_ext(filament_spin_indx_vect,:);
              
              filament_struct_ref(zaehler_f_ref).ccc_filament  = filament_struct(i).ccc_filament_ext(filament_spin_indx_vect,:);
              
              filament_struct_ref(zaehler_f_ref).segment_sensitivity_filament = filament_struct(i).segment_sensitivity_filament_ext(filament_spin_indx_vect,:);
              
              filament_struct_ref(zaehler_f_ref).vect_filament = filament_struct(i).vect_filament_ext(filament_spin_indx_vect,:);
              
              filament_struct_ref(zaehler_f_ref).actin_data_matrix = filament_struct(i).actin_data_matrix_ext(filament_spin_indx_vect,:);
              
              if ~isempty(find(filament_spin_vect(1,filament_spin_indx_vect) == 0.5))
                  error('No valid segment spin!');
              end
              
              filament_struct_ref(zaehler_f_ref).spin_filament = filament_spin_vect(1,filament_spin_indx_vect)';%Does not contain 0.5!
              
              zaehler_f_ref = zaehler_f_ref + 1;
        
        end
        
    end
    
    disp('----------------');
    disp(['Tomogram index: ' num2str(k)]);
    disp(['Number of segments with spin 0: ' num2str(zaehler_seg_tomo_spin_0)]);
    disp(['Percentage of segments with spin 0: ' num2str(round(100.*(zaehler_seg_tomo_spin_0./(zaehler_seg_tomo_spin_0+zaehler_seg_tomo_spin_1))))]);
    disp(['Number of segments with spin 1: ' num2str(zaehler_seg_tomo_spin_1)]);
    disp(['Percentage of segments with spin 1: ' num2str(round(100.*(zaehler_seg_tomo_spin_1./(zaehler_seg_tomo_spin_0+zaehler_seg_tomo_spin_1))))]);
    disp('----------------');
    
    % Reorder filament segments
    for i=1:size(filament_struct_ref,2)
    
        D = pdist2(mean(filament_struct_ref(i).cor_filament),filament_struct_ref(i).cor_filament);
        
        [~,indx_cm] = max(D);
        
        D = pdist2(filament_struct_ref(i).cor_filament(indx_cm,:),filament_struct_ref(i).cor_filament);
        
        [~,indx_sort] = sort(D);
        
        
        filament_struct_ref(i).cor_filament  = filament_struct_ref(i).cor_filament(indx_sort,:);
        
        filament_struct_ref(i).psi_filament  = filament_struct_ref(i).psi_filament(indx_sort,:);
        
        filament_struct_ref(i).ccc_filament  = filament_struct_ref(i).ccc_filament(indx_sort,:);
        
        filament_struct_ref(i).segment_sensitivity_filament = filament_struct_ref(i).segment_sensitivity_filament(indx_sort,:);
        
        filament_struct_ref(i).vect_filament = filament_struct_ref(i).vect_filament(indx_sort,:);
        
        filament_struct_ref(i).actin_data_matrix = filament_struct_ref(i).actin_data_matrix(indx_sort,:);
        
        filament_struct_ref(i).spin_filament = filament_struct_ref(i).spin_filament(indx_sort,:);
        
    end
    
    % Calculate spin properties and confidence values
    for i=1:size(filament_struct_ref,2)
        
        filament_struct_ref(i).spin_filament_mean = mean(filament_struct_ref(i).spin_filament);
        
        if filament_struct_ref(i).spin_filament_mean <= 0.5
            
            filament_struct_ref(i).spin_filament_red = 0;
            filament_struct_ref(i).spin_filament_red_vect = zeros(size(filament_struct_ref(i).cor_filament,1),1);
            
        elseif filament_struct_ref(i).spin_filament_mean > 0.5
            
            filament_struct_ref(i).spin_filament_red = 1;
            filament_struct_ref(i).spin_filament_red_vect = ones(size(filament_struct_ref(i).cor_filament,1),1);
            
        end
        
        % Calculate majority confidence parameter
        filament_struct_ref(i).num_of_seg = size(filament_struct_ref(i).spin_filament,1);
        filament_struct_ref(i).num_of_0 = size(find(filament_struct_ref(i).spin_filament==0),1);
        filament_struct_ref(i).num_of_1 = size(find(filament_struct_ref(i).spin_filament==1),1);
        
        filament_struct_ref(i).frac_of_0 = filament_struct_ref(i).num_of_0./filament_struct_ref(i).num_of_seg;
        filament_struct_ref(i).frac_of_1 = filament_struct_ref(i).num_of_1./filament_struct_ref(i).num_of_seg;
        
        if filament_struct_ref(i).frac_of_0 > filament_struct_ref(i).frac_of_1 %That means filament_struct_ref(i).frac_of_0 > 0.5
            
            % Filament both directions / Majority in spin vector is 0
            filament_struct_ref(i).spin_filament_red_red = 0;
            filament_struct_ref(i).spin_filament_conf = filament_struct_ref(i).frac_of_0;
            
        elseif filament_struct_ref(i).frac_of_0 == filament_struct_ref(i).frac_of_1 %That is 0.5!
            
            % Filament both directions / There is no majority in spin vector
            filament_struct_ref(i).spin_filament_red_red = 0;
            filament_struct_ref(i).spin_filament_conf = filament_struct_ref(i).frac_of_0;
            
        elseif filament_struct_ref(i).frac_of_0 < filament_struct_ref(i).frac_of_1 %That means filament_struct_ref(i).frac_of_1 > 0.5
            
            % Filament both directions / Majority in spin vector is 1
            filament_struct_ref(i).spin_filament_red_red = 1;
            filament_struct_ref(i).spin_filament_conf = filament_struct_ref(i).frac_of_1;
            
        end
        
        % Calculate segment sensitivity confidence
        filament_struct_ref(i).frac_of_pos_sense = size(find(filament_struct_ref(i).segment_sensitivity_filament==1),1)./size(filament_struct_ref(i).spin_filament,1);
        
        % Calculate mean filament ccc
        filament_struct_ref(i).mean_filament_ccc = mean(filament_struct_ref(i).ccc_filament); 
        
        % Calculate combined confidence score / majority and sensitivity
        % filament_struct_ref(i).spin_filament_conf is [0.5,...,1.0]
        % filament_struct_ref(i).frac_of_pos_sense is  [0.0,...,1.0]
        filament_struct_ref(i).filament_conf_com_score = filament_struct_ref(i).spin_filament_conf .* filament_struct_ref(i).frac_of_pos_sense;
        
        % Calculate maximum likelihood estimate of the filaments polarity
        filament_struct_ref(i).filament_mle_polarity = filament_struct_ref(i).spin_filament_conf;
        
        % Calculate standard error of the maximum likelihood estimate
        filament_struct_ref(i).filament_std_error_polarity = sqrt((filament_struct_ref(i).filament_mle_polarity.*(1-filament_struct_ref(i).filament_mle_polarity))./filament_struct_ref(i).num_of_seg);
        
        % Calculate 95% confidence interval for the maximum likelihood estimate
        filament_struct_ref(i).filament_mle_conf_int_a = filament_struct_ref(i).filament_mle_polarity - 1.96.*filament_struct_ref(i).filament_std_error_polarity;
        filament_struct_ref(i).filament_mle_conf_int_b = filament_struct_ref(i).filament_mle_polarity + 1.96.*filament_struct_ref(i).filament_std_error_polarity;
        
    end
    
    % Select filaments based on combined confidence score
    if filament_com_conf_thres ~= -1
        filament_struct_ref_conf = [filament_struct_ref.filament_conf_com_score];
        indx_bad = find(filament_struct_ref_conf < filament_com_conf_thres);
        disp(['Bundle: ' num2str(k) ' | Out of ' num2str(size(filament_struct_ref,2)) ' filaments in this bundle ' num2str(size(indx_bad,2)) ' will be excluded (ccs) | Number of filaments: ' num2str(size(filament_struct_ref,2)-size(indx_bad,2))]);
        filament_struct_ref(indx_bad) = [];
        clear filament_struct_ref_conf;
        clear indx_bad;
    end
    
    % Select filaments based on the 95% confidence interval for the maximum likelihood estimate
    if filament_mle_conf_int_a_thres ~= -1
        filament_struct_ref_filament_mle_conf_int_a = [filament_struct_ref.filament_mle_conf_int_a];
        indx_bad = find(filament_struct_ref_filament_mle_conf_int_a < filament_mle_conf_int_a_thres);
        disp(['Bundle: ' num2str(k) ' | Out of ' num2str(size(filament_struct_ref,2)) ' filaments in this bundle ' num2str(size(indx_bad,2)) ' will be excluded (mle_conf_int) | Number of filaments: ' num2str(size(filament_struct_ref,2)-size(indx_bad,2))]);
        filament_struct_ref(indx_bad) = [];
        clear filament_struct_ref_conf;
        clear indx_bad;
    end
    
    % Extract all filament confidence values
    for i=1:size(filament_struct_ref,2)
    
        spin_filament_conf_all = [spin_filament_conf_all;filament_struct_ref(i).spin_filament_conf];
        
        if filament_struct_ref(i).spin_filament_red_red == 0
             spin_filament_conf_all_0 = [spin_filament_conf_all_0;filament_struct_ref(i).spin_filament_conf];
        elseif filament_struct_ref(i).spin_filament_red_red == 1
             spin_filament_conf_all_1 = [spin_filament_conf_all_1;filament_struct_ref(i).spin_filament_conf];
        end
        
        spin_filament_conf_seg_sense_all = [spin_filament_conf_seg_sense_all;filament_struct_ref(i).frac_of_pos_sense];
        
        spin_filament_conf_mean_ccc_all = [spin_filament_conf_mean_ccc_all;filament_struct_ref(i).mean_filament_ccc];
        
        spin_filament_conf_com_score_all = [spin_filament_conf_com_score_all;filament_struct_ref(i).filament_conf_com_score];
        
        spin_filament_mle_conf_int_a_all = [spin_filament_mle_conf_int_a_all;filament_struct_ref(i).filament_mle_conf_int_a];
        
        spin_filament_conf_tomo_all = [spin_filament_conf_tomo_all;k];
        
    end
    
    % Evaluate segment neighbourhood
    for i=1:size(filament_struct_ref,2)
        
        segment_neighbourhood = zeros(size(filament_struct_ref(i).cor_filament,1),size(filament_struct_ref,2),3);
        
        for m=1:size(filament_struct_ref(i).cor_filament,1)
            
            segment_cor = filament_struct_ref(i).cor_filament(m,:);
            
            for n=1:size(filament_struct_ref,2)
                
                filament_cor = filament_struct_ref(n).cor_filament;
                     
                D = pdist2(segment_cor,filament_cor);
                     
                [val_min_seg,indx_min_seg] = min(D);
                
                segment_neighbourhood(m,n,1) = val_min_seg;
                segment_neighbourhood(m,n,2) = indx_min_seg;
                segment_neighbourhood(m,n,3) = filament_struct_ref(n).spin_filament_red_vect(indx_min_seg,1);% An actin filament has only one direction
                
            end
            
        end
        
        filament_struct_ref(i).segment_neighbourhood = segment_neighbourhood;
        
    end
    
    % Restrict segment neighbourhood to an epsilon-sphere
    for i=1:size(filament_struct_ref,2)
        
        for m=1:size(filament_struct_ref(i).cor_filament,1)
            
            segment_neighbourhood_red = squeeze(filament_struct_ref(i).segment_neighbourhood(m,:,:));
            
            indx_filament = find(segment_neighbourhood_red(:,1) > 0 & segment_neighbourhood_red(:,1) <= epsilon_sphere_radius);
            
            indx_segment = segment_neighbourhood_red(indx_filament,2);
            
            dist_neighbourhood = segment_neighbourhood_red(indx_filament,1);
            
            spin_neighbourhood = segment_neighbourhood_red(indx_filament,3);
            
            filament_struct_ref(i).segment_neighbourhood_epsilon(m).indx_filament = indx_filament';
            filament_struct_ref(i).segment_neighbourhood_epsilon(m).indx_segment = indx_segment';
            filament_struct_ref(i).segment_neighbourhood_epsilon(m).dist_neighbourhood = dist_neighbourhood';
            filament_struct_ref(i).segment_neighbourhood_epsilon(m).spin_neighbourhood = spin_neighbourhood';
            
        end
        
    end
    
    % Calculate filament length and mean distance between segments
    for i=1:size(filament_struct_ref,2)
        
        length_filament_vect = zeros(size(filament_struct_ref(i).cor_filament,1)-1,1);
        length_filament_pre = 0;
        
        for j=1:size(filament_struct_ref(i).cor_filament,1)-1
             
              length_filament_vect(j,1) = pdist2([filament_struct_ref(i).cor_filament(j,  1) filament_struct_ref(i).cor_filament(j,  2) filament_struct_ref(i).cor_filament(j,  3)],...
                                                 [filament_struct_ref(i).cor_filament(j+1,1) filament_struct_ref(i).cor_filament(j+1,2) filament_struct_ref(i).cor_filament(j+1,3)]);
              
              length_filament_pre = length_filament_pre + length_filament_vect(j,1);
            
        end
        filament_struct_ref(i).length_filament_vect = length_filament_vect;
        filament_struct_ref(i).length_filament = length_filament_pre;% multiply by 4 * 3.443 to have the filament length in A
        
    end
    
    % Save refined filaments
    save(FilamentStruct{k},'filament_struct_ref');
    
    %disp(k);
    
end

% Save all bundles distributions
save('/home/Medalia/Projects7/Bruno/ActinPolarity/mapping/mapping_144/bundles/all_bundles_distributions_144.mat',...
    'psi_distribution_all',...
    'ccc_distribution_all',...
    'segment_sensitivity_distribution_all',...
    'tomo_distribution_all',...
    'spin_filament_conf_all',...
    'spin_filament_conf_uni_all',...
    'spin_filament_conf_all_0',...
    'spin_filament_conf_all_1',...
    'spin_filament_conf_seg_sense_all',...
    'spin_filament_conf_mean_ccc_all',...
    'spin_filament_conf_com_score_all',...
    'spin_filament_mle_conf_int_a_all',...
    'spin_filament_conf_tomo_all');

