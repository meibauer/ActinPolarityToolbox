
% Generate BasePath
actin_polarity_construct_base_path;

% Define microscope parameter
Objectpixelsize = 0.3443;%nm

% Check particle coordinates
StartPath = pwd;
particle_avg_check = zeros(36,36,36);
for k=1:size(BasePath,2)
    disp('-----------------------------------------------------------------------');
    disp(['tomogram number ' num2str(k)]);
    cd(BasePath{k});
    overview_wbp = tom_emreadc3f('overview_wbp.em');
    try
    
        load([BasePath{k} '/cor/actin_polarity_cor.mat']);
        
        plist = actin_3d_cor_all;
        plist_filaments = actin_filament_list;
        
        save([BasePath{k} '/cor/actin_polarity_cor_final_144.mat'],'plist','plist_filaments');
        
        disp(['number of particles in this list: ' num2str(size(plist,1))]);
    catch
        disp('contains no particles in this list');
        continue;
    end
    
    del_indx = [];
    zaehler = 1;
    for i=1:size(plist,1)
         try
              particle = overview_wbp(plist(i,1)-18+1:plist(i,1)+18,plist(i,2)-18+1:plist(i,2)+18,plist(i,3)-18+1:plist(i,3)+18);
              particle_avg_check = particle_avg_check + particle;
         catch
              disp(['tomogram ' num2str(k) ' / particle ' num2str(i) ' delete']);
              del_indx(zaehler) = i;
              zaehler = zaehler + 1;
         end
    end
    
    if ~isempty(del_indx)
        
        plist(del_indx,:) = [];
        plist_filaments(del_indx,:) = [];
        
        save([BasePath{k} '/cor/actin_polarity_cor_final_144.mat'],'plist','plist_filaments');
    end
    clear zaehler;
    clear del_indx;
    cd(StartPath);
disp('xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx');
end
clear StartPath i k overview_wbp particle plist;

