
% Construct BasePath
actin_polarity_construct_base_path;

% Create backup folders and move data
for k=1:size(BasePath,2)
    
    %mkdir([BasePath{k} '/segmentation']);
    
    cd([BasePath{k} '/proj/sorted']);
    
    ls
    
    disp(k);
    
    %unix(['mv *.mrc ' BasePath{k} '/segmentation']);
    
end

% % Create backup folders and move data
% for k=1:size(BasePath,2)
%     
%     %mkdir([BasePath{k} '/segmentation']);
%     
%     cd([BasePath{k}]);
%     
%     unix(['mv *.mrc ' BasePath{k} '/segmentation']);
%     
% end

% % Rename previous CTF corrected projections and create folder for new ones
% for k=1:size(BasePath,2)
%     
%     %mkdir([BasePath{k} '/proj/ctfcor_flip_p']);
%     %mkdir([BasePath{k} '/proj/ctfcor_flip_n']);
%     
%     mkdir([BasePath{k} '/proj/prep_flip_p']);
%     mkdir([BasePath{k} '/proj/prep_flip_n']);
%     
%     %unix(['mv ' BasePath{k} '/proj/ctfcor ' BasePath{k} '/proj/ctfcor_mean']);
%     
% end

% % Clear folders
% for k=1:size(BasePath,2)
%     
%     cd([BasePath{k} '/proj/ctfcor_flip_p']);
%     delete *.*;
%     cd([BasePath{k} '/proj/ctfcor_flip_n']);
%     delete *.*;
%     
%     cd([BasePath{k} '/proj/prep_flip_p']);
%     delete *.*;
%     cd([BasePath{k} '/proj/prep_flip_n']);
%     delete *.*;
%     
% end

% % Prepare folders for IMOD reconstruction
% for k=1:size(BasePath,2)
%     
%     mkdir([BasePath{k} '/imod/ctfcor_flip_p']);
%     
% end

% % Add folders for segmentation
% for k=1:size(BasePath,2)
%     
%     %mkdir([BasePath{k} '/imod/segmentation']);
%     mkdir([BasePath{k} '/imod/segmentation/training']);
%     
% end

% % Prepare folders for defocus calibration
% for k=1:size(BasePath,2)
%     
%     mkdir([BasePath{k} '/ps/defocus_calibration/defocus_planes']);
%     mkdir([BasePath{k} '/ps/defocus_calibration/defocus_model_flip_p']);
%     mkdir([BasePath{k} '/ps/defocus_calibration/defocus_model_flip_n']);
%     
% end

