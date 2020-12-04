% Polarity z height analysis

%% Calculate polarity traces

% Define extracted filaments
for t=1:7
    FilamentStruct{t} = ['/home/Medalia/Projects7/Bruno/ActinPolarity/mapping/mapping_144/filaments/filaments_144_tomo_' num2str(t) '.mat'];
end

% Analyze bundles
for k=[1,2,3,4,5,6,7]

    % Load filament structure
    load(FilamentStruct{k},'filament_struct_ref');
    
    % Initialize filament bundle matrix
    filament_bundle_matrix = zeros(size(filament_struct_ref,2),23);
    
    % Calculate center of mass of each filament and extract filament spin
    for i=1:size(filament_struct_ref,2)
         
         filament_bundle_matrix(i,1) = i;
         filament_bundle_matrix(i,2:4) = mean(filament_struct_ref(i).cor_filament,1);
         filament_bundle_matrix(i,5) = filament_struct_ref(i).spin_filament_red_red;
         
    end
    
    % Fit bundle plane
    [X,Y,Z] = actin_polarity_fit_bundle_plane(filament_bundle_matrix(:,2),filament_bundle_matrix(:,3),filament_bundle_matrix(:,4),1024);
    
    % Calculate projection of each filaments center of mass onto the bundle plane
    for i=1:size(filament_struct_ref,2)
         
         d_plane_filament_cm = zeros(1024,1024);
         
         for m=1:1024
              for n=1:1024
              
              d_plane_filament_cm(m,n) = sqrt((filament_bundle_matrix(i,2)-X(m,n)).^2+(filament_bundle_matrix(i,3)-Y(m,n)).^2+(filament_bundle_matrix(i,4)-Z(m,n)).^2); 
              
              end
         end
         
         d_plane_filament_cm_min = min(d_plane_filament_cm(:));
         
         [a,b] = find(d_plane_filament_cm == d_plane_filament_cm_min);
         a = a(1);
         b = b(1);
         
         filament_bundle_matrix(i,6) = a;
         
         filament_bundle_matrix(i,7) = b;
         
         filament_bundle_matrix(i,8) = d_plane_filament_cm(a,b);
         
         filament_bundle_matrix(i,9) = filament_bundle_matrix(i,4)-Z(a,b);
         
    end
    
%     % Plot bundle plane and coordinates
%     
%     figure(k);hold on;plot3(X,Y,Z,'g.');
%     
%     for i=1:size(filament_struct_ref,2)
%     
%          if filament_bundle_matrix(i,9) <= 0
%             
%               plot3(filament_bundle_matrix(i,2),filament_bundle_matrix(i,3),filament_bundle_matrix(i,4),'b+');
%         
%          elseif filament_bundle_matrix(i,9) > 0
%              
%               plot3(filament_bundle_matrix(i,2),filament_bundle_matrix(i,3),filament_bundle_matrix(i,4),'b+');
%         
%          end
%     
%     end
%     
%     view(180,0);xlim([1 1024]);ylim([1 1024]);zlim([1 1024]);
%     xlabel('x-axis');
%     ylabel('y-axis');
%     zlabel('z-axis');
    
    % Find filament center of mass that has greatest distance to the plane
    [~,d_plane_filament_cm_max] = max(filament_bundle_matrix(:,8));
    
    % Create boundary plane by translation plane to that point
    X_boundary_A = X;
    Y_boundary_A = Y;
    Z_boundary_A = Z + sign(filament_bundle_matrix(d_plane_filament_cm_max,9)).*filament_bundle_matrix(d_plane_filament_cm_max,8);
    
    % Find the opposite boundary plane
    
    % Calculate projection of each filaments center of mass onto the bundle plane
    for i=1:size(filament_struct_ref,2)
        
         d_plane_filament_cm = zeros(1024,1024);
         
         for m=1:1024
              for n=1:1024
                
                   d_plane_filament_cm(m,n) = sqrt((filament_bundle_matrix(i,2)-X_boundary_A(m,n)).^2+(filament_bundle_matrix(i,3)-Y_boundary_A(m,n)).^2+(filament_bundle_matrix(i,4)-Z_boundary_A(m,n)).^2);
                
              end
         end
         
         d_plane_filament_cm_min = min(d_plane_filament_cm(:));
         
         [a,b] = find(d_plane_filament_cm == d_plane_filament_cm_min);
         a = a(1);
         b = b(1);
         
         filament_bundle_matrix(i,10) = a;
         
         filament_bundle_matrix(i,11) = b;
         
         filament_bundle_matrix(i,12) = d_plane_filament_cm(a,b);
         
         filament_bundle_matrix(i,13) = filament_bundle_matrix(i,4)-Z_boundary_A(a,b);
         
    end
    
    % Find filament center of mass that has the greatest distance to that plane
    [~,d_plane_filament_cm_max] = max(filament_bundle_matrix(:,12));
    
%     % Create boundary plane by translation plane to that point
%     X_boundary_B = X_boundary_A;
%     Y_boundary_B = Y_boundary_A;
%     Z_boundary_B = Z_boundary_A + sign(filament_bundle_matrix(d_plane_filament_cm_max,13)).*filament_bundle_matrix(d_plane_filament_cm_max,12);
    
    % Calculate maximum plane translation that is the distance between the two boundary planes
    plane_translation_max = round(sign(filament_bundle_matrix(d_plane_filament_cm_max,13)).*filament_bundle_matrix(d_plane_filament_cm_max,12));
    
    % Initialize z height analysis
    z_height_stat = zeros(size(0:sign(plane_translation_max).*8:plane_translation_max,2),6);
    
    % Peform z height analysis
    zaehler = 1;
    for p_step=0:sign(plane_translation_max).*8:plane_translation_max
         
         % Translate moving plane to new position
         X_plane_new = X_boundary_A;
         Y_plane_new = Y_boundary_A;
         Z_plane_new = Z_boundary_A + p_step;
         
         % Calculate projection of each filaments center of mass onto the bundle plane
         for i=1:size(filament_struct_ref,2)
            
              d_plane_filament_cm = zeros(1024,1024);
              
              for m=1:1024
                   for n=1:1024
                    
                        d_plane_filament_cm(m,n) = sqrt((filament_bundle_matrix(i,2)-X_plane_new(m,n)).^2+(filament_bundle_matrix(i,3)-Y_plane_new(m,n)).^2+(filament_bundle_matrix(i,4)-Z_plane_new(m,n)).^2);
                    
                   end
              end
              
              d_plane_filament_cm_min = min(d_plane_filament_cm(:));
              
              [a,b] = find(d_plane_filament_cm == d_plane_filament_cm_min);
              a = a(1);
              b = b(1);
              
              filament_bundle_matrix(i,20) = a;
              
              filament_bundle_matrix(i,21) = b;
              
              filament_bundle_matrix(i,22) = d_plane_filament_cm(a,b);
              
              filament_bundle_matrix(i,23) = filament_bundle_matrix(i,4)-Z_plane_new(m,n);
            
         end
         
         % Find candidate points to current plane
         indx_points = find(filament_bundle_matrix(:,22)<=12);
         
         % Extract spins
         spin_points = filament_bundle_matrix(indx_points,5);
         
         z_height_stat(zaehler,1) = size(find(spin_points==0),1);%spin=0
         z_height_stat(zaehler,2) = size(find(spin_points==1),1);%spin=1
         
         z_height_stat(zaehler,3) = z_height_stat(zaehler,1) ./ (z_height_stat(zaehler,1) + z_height_stat(zaehler,2));
         z_height_stat(zaehler,4) = z_height_stat(zaehler,2) ./ (z_height_stat(zaehler,1) + z_height_stat(zaehler,2));
         
         z_height_stat(zaehler,5) = p_step;
         z_height_stat(zaehler,6) = mean(Z_plane_new(:));
         
         zaehler = zaehler + 1;
         
    end
    
    % Store z height analysis in structure
    actin_polarity_layers(k).filament_bundle_matrix = filament_bundle_matrix;
    actin_polarity_layers(k).z_height_stat = z_height_stat;
    
    disp(k);
    
end

% Save z height analysis
save('/home/Medalia/Projects7/Bruno/ActinPolarity/mapping/mapping_144/z_height_analysis/workspace_z_height_analysis.mat','actin_polarity_layers');
clear all;


%% Analyze polarity traces

% Load z height analysis structure
load('/home/Medalia/Projects7/Bruno/ActinPolarity/mapping/mapping_144/z_height_analysis/workspace_z_height_analysis.mat','actin_polarity_layers');

% Define bundles that need to be flipped
indx_bundles_flip = [4 5];

% Correct polarity traces
for k=indx_bundles_flip
    
    actin_polarity_layers(k).z_height_stat = actin_polarity_layers(k).z_height_stat(:,[2 1 4 3 5 6]);
    
end

% Remove planes with less than 6 filaments
% Sort starting from support
% Normalize bundle thickness
for k=1:7
    
    % Extract z height statistics
    z_height_stat = actin_polarity_layers(k).z_height_stat;
    
    if ~isempty(z_height_stat)
        
        % Remove planes with less than 6 filaments
        filaments_sum = sum(z_height_stat(:,1:2),2);
        indx_del_z_height_stat = find(sum(z_height_stat(:,1:2),2) < 6);
        
        % Delete those planes
        z_height_stat(indx_del_z_height_stat,:) = [];
        
        % Sort starting from support
        [~,indx_sort_support] = sort(z_height_stat(:,6));
        z_height_stat_prep = z_height_stat(indx_sort_support,:);
        
        % Extract estimated thickness of bundle
        z_dim_bundle = max(z_height_stat_prep(:,6))-min(z_height_stat_prep(:,6));
        z_dim_bundle_nm = (4.*3.443/10).*(max(z_height_stat_prep(:,6))-min(z_height_stat_prep(:,6)));
        
        % Normalize bundle thickness
        z_height_stat_prep(:,7) = tom_norm(z_height_stat_prep(:,6),1);
        
        % Interpolate polarity trace / spin=0
        z_polarity_trace_spin_0 = interp1(z_height_stat_prep(:,7),z_height_stat_prep(:,3),0:0.1:1,'spline');
        z_polarity_trace_spin_1 = interp1(z_height_stat_prep(:,7),z_height_stat_prep(:,4),0:0.1:1,'spline');
        
        % Calculate difference
        z_polarity_trace_diff = z_polarity_trace_spin_0 - z_polarity_trace_spin_1;
        
        % Place back in z height statistics structure
        actin_polarity_layers(k).z_height_stat = z_height_stat_prep;
        actin_polarity_layers(k).z_dim_bundle = z_dim_bundle;
        actin_polarity_layers(k).z_dim_bundle_nm = z_dim_bundle_nm;
        
        actin_polarity_layers(k).z_polarity_trace_x = 0:0.1:1;
        actin_polarity_layers(k).z_polarity_trace_spin_0 = z_polarity_trace_spin_0;
        actin_polarity_layers(k).z_polarity_trace_spin_1 = z_polarity_trace_spin_1;
        
        actin_polarity_layers(k).z_polarity_trace_diff = z_polarity_trace_diff;
        
    end
    
end

% Plot interpolated actin polarity traces

figure(1);hold on;

zaehler = 1;

for k=1:7
    
    if ~isempty(actin_polarity_layers(k).z_height_stat)
    
         subplot(1,7,zaehler);hold on;plot(actin_polarity_layers(k).z_polarity_trace_spin_0,actin_polarity_layers(k).z_polarity_trace_x,'b-');
         subplot(1,7,zaehler);hold on;plot(actin_polarity_layers(k).z_polarity_trace_spin_1,actin_polarity_layers(k).z_polarity_trace_x,'r-');
         box on;
         
         xlim([0 1]);xticks([0 0.5 1]);
         ylim([0 1]);xticks([0 0.5 1]);
         
         title(k);
         xlabel('Spin ratio');
         ylabel('Position');
         
         zaehler = zaehler + 1;
         
    end
    
end

% % Plot interpolated actin polarity trace differences
% 
% figure(2);hold on;
% 
% zaehler = 1;
% 
% for k=1:7
%     
%     if ~isempty(actin_polarity_layers(k).z_height_stat)
%     
%          subplot(1,7,zaehler);hold on;plot(actin_polarity_layers(k).z_polarity_trace_x,actin_polarity_layers(k).z_polarity_trace_diff,'g-');
%          subplot(1,7,zaehler);hold on;plot(actin_polarity_layers(k).z_polarity_trace_x,zeros(1,size(actin_polarity_layers(k).z_polarity_trace_x,2)),'k--');
%          
%          xlim([0 1]);xticks([0 0.5 1]);
%          ylim([-1 +1]);xticks([0 0.5 1]);
%          
%          title(k);
%          xlabel('Position');
%          ylabel('Spin ratio difference');
%          
%          zaehler = zaehler + 1;
%          
%     end
%     
% end

