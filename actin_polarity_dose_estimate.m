
% Generate BasePath
actin_polarity_construct_base_path;

% Define microscope parameter
Objectpixelsize = 0.3443;%nm
dose_at_zero_deg = 1.2;%e-/A^2

% Estimate cumulative dosage
for k=1:1%size(BasePath,2)
    
    % Load tiltangles
    load([BasePath{k} '/alg/alg.mat'],'alg0');
    tiltangles = alg0.Tiltangles;
    
    % Estimate dose per tilt
    dose_per_tilt = dose_at_zero_deg .* (1./cosd(abs(tiltangles)));%e-/A^2
    
    % Initialize cumulative dose
    cumulative_dose = zeros(1,size(tiltangles,2));
    
    % Find start projection
    [~,ind_start_proj] = min(abs(tiltangles + 30));
    
    % Initialize the sum
    sum_dose = 0;
    
    for i=ind_start_proj:+1:size(tiltangles,2)
    
        sum_dose = dose_per_tilt(1,i) + sum_dose;
        cumulative_dose(1,i) = sum_dose;

    end
    
    for i=(ind_start_proj-1):-1:1
    
        sum_dose = dose_per_tilt(1,i) + sum_dose;
        cumulative_dose(1,i) = sum_dose;
    
    end
    
end

% Plot dose estimate curve
figure(1);hold on;
plot(cumulative_dose,'b-','LineWidth',1);
plot(40.*ones(1,61),'k--');
box off;
xlim([1 61]);ylim([0 100]);
set(gca,'XTick',[1 16  31 46 61]);
set(gca,'XTickLabel',{'-60' '-30' '0' '+30' '+60'});
xlabel('Projection index');
ylabel('Cumulative electron dose (e-/A^2)');

