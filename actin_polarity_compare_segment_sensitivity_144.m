
% Load confidence distribution and tomo membership of all segments
load('/home/Medalia/Projects7/Bruno/ActinPolarity/mapping/mapping_144/bundles/all_bundles_distributions_144.mat','spin_filament_conf_tomo_all','spin_filament_conf_seg_sense_all');

% % Plot all psi histograms
% for k=1:7
%    
%     indx_sel = find(tomo_distribution_all == k);
%     
%     psi_sel = psi_distribution_all(indx_sel);
%     
%     figure(1);subplot(1,7,k);histogram(psi_sel,360);xlim([-200 +200]);title(['Psi (' num2str(k) ')']);
%     
% end

% Estimate bundle direction
bundle_dir = abs([-38.1     -1.0     +55.9     +71.2     +70.3     -67.8     -49.2]);

%figure(2);histogram(bundle_dir,10);xlim([0 90]);

% Create all bundle direction
bundle_dir_all = zeros(size(spin_filament_conf_tomo_all,1),1);
for k=1:7
   
    indx_sel = find(spin_filament_conf_tomo_all == k);
    
    bundle_dir_all(indx_sel) = bundle_dir(k);
    
end

figure;figure_indx = gcf;figure(figure_indx);plot(bundle_dir_all,spin_filament_conf_seg_sense_all,'r+');

% Calculate mean segment sensitivity confidence per bundle direction
bundle_dir_conf = zeros(1,7);
for k=1:7
   
    indx_sel = find(spin_filament_conf_tomo_all == k);
    
    bundle_dir_conf(1,k) = mean(spin_filament_conf_seg_sense_all(indx_sel));
    bundle_dir_conf_std(1,k) = std(spin_filament_conf_seg_sense_all(indx_sel));
    
end

figure(figure_indx);hold on;plot(bundle_dir,bundle_dir_conf,'b+','MarkerSize',20,'LineWidth',2);

xlim([-5 77]);
ylim([0.45 1.05]);
xlabel('Angle tiltaxis with bundle orientation (deg)');
ylabel('Sensitivity confidence score');

% Fit mean filament conf per bundle values
P = robustfit(rmmissing(bundle_dir),rmmissing(bundle_dir_conf));
yfit = P(1)+P(2).*rmmissing(bundle_dir);

figure(figure_indx);plot(rmmissing(bundle_dir),yfit,'b-','LineWidth',2);

% Add 0.80 (=16/20) threshold to figure
figure(figure_indx);plot(-5:1:90,0.80.*ones(1,96),'k--');

% Plot histogram of confidence score values
figure;figure_indx = gcf;figure(figure_indx);histogram(spin_filament_conf_seg_sense_all,'BinEdges',[0.0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0]);
xlabel('Sensitivity confidence score');
ylabel('Number of rec. filaments');
ylim([0 1900]);

% Number and percentage of filaments above threshold
size(find(spin_filament_conf_seg_sense_all>=0.80),1)
size(find(spin_filament_conf_seg_sense_all>=0.80),1)./size(find(spin_filament_conf_seg_sense_all>=0.0),1)
% Length threshold = 3 ---> ~78% (2252 of 2893 filaments) of all filaments with a length of at least 3 segments are equal or better this score

