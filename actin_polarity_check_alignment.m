
% Generate BasePath
actin_polarity_construct_base_path;

% Extract resval values
for k=1:size(BasePath,2)

    % Load alignment
    load([BasePath{k} '/alg/alg.mat'],'alg0','algXY');
    
    % Extract original tiltangles
    tiltangles_raw = alg0.mf(1,:,1);
    
    % Extract tiltangles and tiltaxis from initial alignment
    tiltangles_init = alg0.Tiltangles;
    tiltaxis_init = alg0.Tiltaxis;
    
    % Extract optimized tiltangles and tiltaxis
    tiltangles_alg = algXY.Tiltangles;
    tiltaxis_alg = algXY.Tiltaxis;
    
    % Calculate statistics and differences
    disp('------------------------------------------------');
    disp(BasePath{k});
    disp(k);
    disp('------------------------------------------------');
    tom_dev(abs(tiltangles_raw - tiltangles_init));%--->identical
    all_tiltaxis_init(1,k) = tom_dev(tiltaxis_init);%---> std ~ 0.15 deg, mean ~ 264 deg
    tiltangles_check(1,k) = tom_dev(abs(tiltangles_raw - tiltangles_alg));%---> only tomogram k=2,6 shows a deviation >0.5 deg
    tom_dev(abs(tiltaxis_init - tiltaxis_alg));%---> small deviation
    tom_dev(tiltaxis_alg);
    disp('------------------------------------------------');
    disp(' ');
    disp(' ');
    disp(' ');
    
    all_resval(1,k) = algXY.resval.*0.3443.*10;%A
    all_marker(1,k) = size(algXY.mf,3);
    
end

figure(1);plot(all_resval,'LineWidth',1);hold on;plot(mean(all_resval).*ones(1,18),'k--');
xlim([1 7]);ylim([2 8]);box off;xlabel('Tomogram index');ylabel('Alignment residual (A)');

figure(2);plot(all_marker,'LineWidth',1);hold on;plot(mean(all_marker).*ones(1,18),'k--');
xlim([1 7]);ylim([0 12]);box off;xlabel('Tomogram index');ylabel('Number of marker');

