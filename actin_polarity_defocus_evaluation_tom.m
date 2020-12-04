
% Generate BasePath
actin_polarity_construct_base_path;

% Define microscope parameter
Objectpixelsize = 0.3443;%nm
Voltage = 300;%kV
Cs = 2.7;%mm

% % Calculate mean power spectrum
% for k=1:size(BasePath,2)
%     PerformMeanPowerSpectrum(BasePath{k},256,2);
% end

% % Evaluate mean power spectra
% for k=1:size(BasePath,2)
%     disp('--------------------------------------------------------------');
%     disp(['Base path: ' BasePath{k} '  / Tomogram index: ' num2str(k)]);
%     disp('--------------------------------------------------------------');
%     load([BasePath{k} '/ps/mean_power_spectrum.mat'],'ps');
%     disp('mean_power_spectrum under consideration ...');
%     [psf_256_2,dz_256_2] = InspectDZgui(ps,Objectpixelsize,Voltage,Cs,-4,-4,0,1,8,0);
%     dz_est(1,k) = dz_256_2;
%     save([BasePath{k} '/ps/mean_defocus_256_2.mat'],'psf_256_2','dz_256_2');
%     clear ps;
% end

% % Calculate mean defocus
% dz_est_mean = mean(dz_est);
% disp(dz_est_mean);%---> -4.64 mum

% % Plot mean defocus values
% figure(1);hold on;
% plot(dz_est,'LineWidth',1);
% plot(mean(dz_est) .* ones(1,7),'k--');
% xlim([1 7]);
% ylim([-6 -4]);
% xlabel('Tomogram index');
% ylabel('Mean defocus (mum)');

% Create CTF parameter structure
dz_in.Objectpixelsize = Objectpixelsize;%nm /// the pixelsize
dz_in.Voltage = Voltage;%kV /// the voltage
dz_in.Cs = Cs;%mum /// the spherical aberration constant / TITAN
dz_in.W = 0;%unity /// the amplitude contrast fraction
dz_in.dz_est = -4.64;%mum /// the estimated mean defocus (eg. autofocus) that controls masking for mean tiltseries defocus fit
dz_in.dz_start = -12;%mum /// the start value of the defocus search for both mean defocus search and per projection defocus search
dz_in.dz_incr = 0.05;%mum /// the increment of the defocus search
dz_in.dz_end = -0.5;%mum /// the end value of the defocus search (dz_in.dz_start < dz_in.dz_end) for both mean defocus search and per projection defocus search
dz_in.dz_scan_interval_astig = 5;%mum /// defocus fit scan range astigmatism
dz_in.fs = 256;%pixel /// fieldsize for periodogram averaging
dz_in.ov = 2;%unity /// overlap factor
dz_in.highp = 1;%pixel /// high-pass power spectra enhancement filter
dz_in.lowp = 8;%pixel /// low-pass power spectra enhancement filter

% % Analyze tiltseries defocus
% for k=1:size(BasePath,2)
%     
%     % Load start value
%     load([BasePath{k} '/ps/mean_defocus_256_2.mat'],'dz_256_2');
%     dz_in.dz_est = dz_256_2;
%     
%     % Analyze mean defocus and mean defocus of each projection
%     [dz_eval_mean,dz_eval_proj] = AnalyzeTiltseriesDefocus([BasePath{k} '/proj/sorted'],'sorted_',[BasePath{k} '/alg'],'alg.mat',dz_in,0,1,7,0);
%     
%     % Save results
%     save([BasePath{k} '/ps/dz_eval_tom.mat'],'dz_eval_mean','dz_eval_proj');
%     
%     disp(k);
%     
% end

% Remove defocus measurement outliers and create defocus correction values
for k=1:size(BasePath,2)
   
    % Load defocus evaluation / tom
    load([BasePath{k} '/ps/dz_eval_tom.mat'],'dz_eval_proj');
    
    % Extract defocus measurements
    dz_eval_proj = [dz_eval_proj.estdz_proj];%tom
    
    % Mean and std of defocus values within search range
    mean_dz_eval_proj = mean(dz_eval_proj(dz_eval_proj<=-1));
    std_dz_eval_proj =   std(dz_eval_proj(dz_eval_proj<=-1));
    
    % Identify outlier defocus values
    indx_bad = find(dz_eval_proj < mean_dz_eval_proj-3.*std_dz_eval_proj | dz_eval_proj > mean_dz_eval_proj+3.*std_dz_eval_proj);
    
    % Identify acceptable defocus values
    indx_good = find(dz_eval_proj >= mean_dz_eval_proj-3.*std_dz_eval_proj & dz_eval_proj <= mean_dz_eval_proj+3.*std_dz_eval_proj);
    
    % Recalculate mean defocus
    mean_dz_eval_proj_mod = mean(dz_eval_proj(indx_good));
    
    % Replace outlier defocus values by mean defocus
    dz_eval_proj_mod = dz_eval_proj;
    dz_eval_proj_mod(indx_bad) = mean_dz_eval_proj_mod;
    clear dz_eval_proj;
    
    % Save defocus correction values
    save([BasePath{k} '/ps/dz_eval_tom_mod.mat'],'dz_eval_proj_mod');
    
    % Plot defocus correction values
    figure(1);
    hold on;
    plot(dz_eval_proj_mod,'b-');
    plot(indx_good,dz_eval_proj_mod(indx_good),'b.');
    plot(indx_bad, dz_eval_proj_mod(indx_bad),'r+');
    plot(mean(dz_eval_proj_mod).*ones(1,size(dz_eval_proj_mod,2)),'b--');
    
    xlim([1 size(dz_eval_proj_mod,2)]);
    ylim([-7 -3]);
    title(['Defocus measurement / tiltseries: ' num2str(k)]);
    box on;
    hold off;
    figure(1);drawnow;pause(0.5);
    saveas(gcf,['/home/Medalia/Projects7/Bruno/ActinPolarity/doc/vis_defocus/actin_polarity_tomo_def_' num2str(k) '.png'],'png');
    delete(gcf);
    
    disp(k);
    
end

