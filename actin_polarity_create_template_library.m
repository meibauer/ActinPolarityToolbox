% Actin polarity mapping / template library

%% Template library / Box size 144 pixels

% Load bundles structure / Refine3D / job017
actin_ref = tom_mrcread('/home/Medalia/Projects8/ActinBundles/actin_bundles_structure_rev/prealignment_3/relion/Refine3D/job017/run_class001.mrc');
actin_ref = actin_ref.Value;

% Rotate structure parallel to x-axis
actin_ref_rot = tom_rotate(actin_ref,[270 90 90]);

% Initialize projections
actin_ref = zeros(144,144,120);

% Projection loop
zaehler = 1;
for k=0:3:360-3
    
    % Rotate structure
    actin_ref_pre = tom_rotate(actin_ref_rot,[0 0 k]);
    
    % Project structure
    actin_ref_pre = mean(actin_ref_pre,3);
    
    % Normalize structure
    actin_ref_pre = (actin_ref_pre - mean(actin_ref_pre(:)))./std(actin_ref_pre(:));
    
    actin_ref(:,:,zaehler) = actin_ref_pre;
    
    disp(zaehler);
    
    zaehler = zaehler + 1;
    
end

% Save reference stack
tom_mrcwrite(single(actin_ref),'name','/home/Medalia/Projects7/Bruno/ActinPolarity/template_library/template_library_prototype/Template_Lib/actin_template_library_box_size_144.mrcs');


%% Template library / Box size 104 pixels

% Load template library
template_lib = tom_mrcreadf('/home/Medalia/Projects7/Bruno/ActinPolarity/template_library/template_library_prototype/Template_Lib/actin_template_library_box_size_144.mrcs');

% Initialize the new template library
template_lib_new = zeros(104,104,120);

% Create the new template library
for k=1:120
   
    % Cut out structure
    template_lib_new_pre = template_lib(72-52+1:72+52,72-52+1:72+52,k);
    
    % Renormalize structure
    template_lib_new_pre = (template_lib_new_pre - mean(template_lib_new_pre(:)))./std(template_lib_new_pre(:));
    
    % Paste into stack
    template_lib_new(:,:,k) = template_lib_new_pre;
    
    disp(k);
    
end

% Save reference stack
tom_mrcwrite(single(template_lib_new),'name','/home/Medalia/Projects7/Bruno/ActinPolarity/template_library/template_library_prototype/Template_Lib/actin_template_library_box_size_104.mrcs');


%% Check templates

% Load the original structure
actin_ref = tom_mrcreadf('/home/Medalia/Projects8/ActinBundles/actin_bundles_structure_rev/prealignment_3/relion/Refine3D/job017/run_class001.mrc');

% Load the prototype template
actin_template_prototype = tom_mrcreadf('/home/Medalia/Projects7/Bruno/ActinPolarity/template_library/template_library_prototype/Template_Lib/run_class001_bundles_job017.mrc');

% Check difference
tom_dev(actin_ref-actin_template_prototype);

% ---> These files are identical

% Load the original structure used for confirmation of cross resolution
actin_ref_cross_res = tom_mrcreadf('/home/Medalia/Projects8/ActinBundles/doc/actin_resolution_check/run_class001_bundles_job017.mrc');

% Check difference
tom_dev(actin_ref-actin_ref_cross_res);

% ---> These files are identical
% ---> The file Template_Lib/run_class001_bundles_job017.mrc is the original left handed Actin structure


%% Add 2d solvent mask for psi fine alignment

% Load original prototype mask
mask_2d = tom_mrcreadf('/home/Medalia/Projects7/Bruno/Actin_MEFs/Manual_Segmentation/relion_proj/mask.mrc');

% Save mask / Box size 144 pixels
tom_mrcwrite(mask_2d,'name','/home/Medalia/Projects7/Bruno/ActinPolarity/template_library/template_library_prototype/Template_Lib/mask_2d_144.mrc');

% Cut mask / Box size 104 pixels
mask_2d_cut = mask_2d(72-52+1:72+52,72-52+1:72+52);

% Save mask / Box size 104 pixels
tom_mrcwrite(mask_2d_cut,'name','/home/Medalia/Projects7/Bruno/ActinPolarity/template_library/template_library_prototype/Template_Lib/mask_2d_104.mrc');


%% Create right handed helical reconstruction templates with orthogonal directions from the original bundles structure / Box size 144 pixels

% Load bundles structure / Refine3D / job017
actin_ref = tom_mrcread('/home/Medalia/Projects8/ActinBundles/actin_bundles_structure_rev/prealignment_3/relion/Refine3D/job017/run_class001.mrc');
actin_ref = actin_ref.Value;

% Save bundles structure
tom_mrcwrite(actin_ref,'name','/home/Medalia/Projects7/Bruno/ActinPolarity/template_library/run_class001_bundles_job017_144.mrc');

% Save bundles structure with reversed direction
tom_mrcwrite(tom_rotate(actin_ref,[0 0 180]),'name','/home/Medalia/Projects7/Bruno/ActinPolarity/template_library/run_class001_bundles_job017_144_rot180.mrc');

% Mirror bundles structure
relion_image_handler --i /home/Medalia/Projects7/Bruno/ActinPolarity/template_library/run_class001_bundles_job017_144.mrc --o /home/Medalia/Projects7/Bruno/ActinPolarity/template_library/run_class001_bundles_job017_144_mirrored.mrc --angpix 3.443 --flipZ

% Mirror bundles structure with reversed direction
relion_image_handler --i /home/Medalia/Projects7/Bruno/ActinPolarity/template_library/run_class001_bundles_job017_144_rot180.mrc --o /home/Medalia/Projects7/Bruno/ActinPolarity/template_library/run_class001_bundles_job017_144_rot180_mirrored.mrc --angpix 3.443 --flipZ

% Apply helical parameter to check whether they are compatible
relion_helix_toolbox --impose --i /home/Medalia/Projects7/Bruno/ActinPolarity/template_library/run_class001_bundles_job017_144_mirrored.mrc --o /home/Medalia/Projects7/Bruno/ActinPolarity/template_library/run_class001_bundles_job017_144_mirrored_sym.mrc --angpix 3.443 --rise 27.6 --twist -166.7 --cyl_outer_diameter 110 --nr_asu 18 --z_percentage 0.3

% Apply helical parameter to check whether they are compatible
relion_helix_toolbox --impose --i /home/Medalia/Projects7/Bruno/ActinPolarity/template_library/run_class001_bundles_job017_144_rot180_mirrored.mrc --o /home/Medalia/Projects7/Bruno/ActinPolarity/template_library/run_class001_bundles_job017_144_rot180_mirrored_sym.mrc --angpix 3.443 --rise 27.6 --twist -166.7 --cyl_outer_diameter 110 --nr_asu 18 --z_percentage 0.3


%% Create right handed helical reconstruction templates with orthogonal directions from the original bundles structure / Box size 104 pixels

% Load bundles structure / Refine3D / job017
actin_ref_pre = tom_mrcread('/home/Medalia/Projects8/ActinBundles/actin_bundles_structure_rev/prealignment_3/relion/Refine3D/job017/run_class001.mrc');
actin_ref_pre = actin_ref_pre.Value;

% Cut box
actin_ref = actin_ref_pre(72-52+1:72+52,72-52+1:72+52,72-52+1:72+52);

% Save bundles structure
tom_mrcwrite(actin_ref,'name','/home/Medalia/Projects7/Bruno/ActinPolarity/template_library/run_class001_bundles_job017_104.mrc');

% Save bundles structure with reversed direction
tom_mrcwrite(tom_rotate(actin_ref,[0 0 180]),'name','/home/Medalia/Projects7/Bruno/ActinPolarity/template_library/run_class001_bundles_job017_104_rot180.mrc');

% Mirror bundles structure
relion_image_handler --i /home/Medalia/Projects7/Bruno/ActinPolarity/template_library/run_class001_bundles_job017_104.mrc --o /home/Medalia/Projects7/Bruno/ActinPolarity/template_library/run_class001_bundles_job017_104_mirrored.mrc --angpix 3.443 --flipZ

% Mirror bundles structure with reversed direction
relion_image_handler --i /home/Medalia/Projects7/Bruno/ActinPolarity/template_library/run_class001_bundles_job017_104_rot180.mrc --o /home/Medalia/Projects7/Bruno/ActinPolarity/template_library/run_class001_bundles_job017_104_rot180_mirrored.mrc --angpix 3.443 --flipZ

% Apply helical parameter to check whether they are compatible
relion_helix_toolbox --impose --i /home/Medalia/Projects7/Bruno/ActinPolarity/template_library/run_class001_bundles_job017_104_mirrored.mrc --o /home/Medalia/Projects7/Bruno/ActinPolarity/template_library/run_class001_bundles_job017_104_mirrored_sym.mrc --angpix 3.443 --rise 27.6 --twist -166.7 --cyl_outer_diameter 110 --nr_asu 13 --z_percentage 0.3

% Apply helical parameter to check whether they are compatible
relion_helix_toolbox --impose --i /home/Medalia/Projects7/Bruno/ActinPolarity/template_library/run_class001_bundles_job017_104_rot180_mirrored.mrc --o /home/Medalia/Projects7/Bruno/ActinPolarity/template_library/run_class001_bundles_job017_104_rot180_mirrored_sym.mrc --angpix 3.443 --rise 27.6 --twist -166.7 --cyl_outer_diameter 110 --nr_asu 13 --z_percentage 0.3


%% Prepare cubic reference map from EMD-6179

% Read reference structure (has a pixel size of 1.05 A)
ref_pre = tom_mrcread('/home/Medalia/Projects7/Bruno/ActinPolarity/template_library/EMD-6179/map/emd_6179.map');
ref_pre = single(ref_pre.Value);

% Put it in cubic box
ref = zeros(200,200,200);
ref(100-80+1:100+80,100-80+1:100+80,:) = ref_pre;

% Save reference in cubic box
tom_mrcwrite(ref,'name','/home/Medalia/Projects7/Bruno/ActinPolarity/template_library/emd_6179_200_pixels.mrc');clear ref_pre;clear ref;


%% Prepare reference map with 6.886 A resolution / Box size 144 pixels

% Rescale reference structure
relion_image_handler --i /home/Medalia/Projects7/Bruno/ActinPolarity/template_library/emd_6179_200_pixels.mrc --o /home/Medalia/Projects7/Bruno/ActinPolarity/template_library/emd_6179_144_pixels_angpix_3_4A_res_6_886A.mrc --angpix 1.05 --rescale_angpix 3.443 --new_box 144 --lowpass 6.886

% Add helical symmetry
relion_helix_toolbox --impose --i /home/Medalia/Projects7/Bruno/ActinPolarity/template_library/emd_6179_144_pixels_angpix_3_4A_res_6_886A.mrc --o /home/Medalia/Projects7/Bruno/ActinPolarity/template_library/emd_6179_144_pixels_angpix_3_4A_res_6_886A_sym.mrc --angpix 3.443 --rise 27.6 --twist -166.7 --cyl_outer_diameter 200 --nr_asu 18 --z_percentage 0.3


%% Prepare reference map with 8 A resolution  / Box size 144 pixels

% Rescale reference structure
relion_image_handler --i /home/Medalia/Projects7/Bruno/ActinPolarity/template_library/emd_6179_200_pixels.mrc --o /home/Medalia/Projects7/Bruno/ActinPolarity/template_library/emd_6179_144_pixels_angpix_3_4A_res_8A.mrc --angpix 1.05 --rescale_angpix 3.443 --new_box 144 --lowpass 8

% Add helical symmetry
relion_helix_toolbox --impose --i /home/Medalia/Projects7/Bruno/ActinPolarity/template_library/emd_6179_144_pixels_angpix_3_4A_res_8A.mrc --o /home/Medalia/Projects7/Bruno/ActinPolarity/template_library/emd_6179_144_pixels_angpix_3_4A_res_8A_sym.mrc --angpix 3.443 --rise 27.6 --twist -166.7 --cyl_outer_diameter 200 --nr_asu 18 --z_percentage 0.3


%% Filter template of rln_144/Refine3D/job006 to 30 A resolution

% Apply lowpass filter
relion_image_handler --i /home/Medalia/Projects7/Bruno/ActinPolarity/template_library/run_class001_bundles_job017_144_mirrored_sym.mrc --o /home/Medalia/Projects7/Bruno/ActinPolarity/template_library/run_class001_bundles_job017_144_mirrored_sym_lowpass_30A.mrc --angpix 3.443 --lowpass 30


%% Prepare reference map with 6.886 A resolution / Box size 104 pixels

% Rescale reference structure
relion_image_handler --i /home/Medalia/Projects7/Bruno/ActinPolarity/template_library/emd_6179_200_pixels.mrc --o /home/Medalia/Projects7/Bruno/ActinPolarity/template_library/emd_6179_104_pixels_angpix_3_4A_res_6_886A.mrc --angpix 1.05 --rescale_angpix 3.443 --new_box 104 --lowpass 6.886

% Add helical symmetry
relion_helix_toolbox --impose --i /home/Medalia/Projects7/Bruno/ActinPolarity/template_library/emd_6179_104_pixels_angpix_3_4A_res_6_886A.mrc --o /home/Medalia/Projects7/Bruno/ActinPolarity/template_library/emd_6179_104_pixels_angpix_3_4A_res_6_886A_sym.mrc --angpix 3.443 --rise 27.6 --twist -166.7 --cyl_outer_diameter 200 --nr_asu 13 --z_percentage 0.3


%% Prepare reference map with 8 A resolution  / Box size 104 pixels

% Rescale reference structure
relion_image_handler --i /home/Medalia/Projects7/Bruno/ActinPolarity/template_library/emd_6179_200_pixels.mrc --o /home/Medalia/Projects7/Bruno/ActinPolarity/template_library/emd_6179_104_pixels_angpix_3_4A_res_8A.mrc --angpix 1.05 --rescale_angpix 3.443 --new_box 104 --lowpass 8

% Add helical symmetry
relion_helix_toolbox --impose --i /home/Medalia/Projects7/Bruno/ActinPolarity/template_library/emd_6179_104_pixels_angpix_3_4A_res_8A.mrc --o /home/Medalia/Projects7/Bruno/ActinPolarity/template_library/emd_6179_104_pixels_angpix_3_4A_res_8A_sym.mrc --angpix 3.443 --rise 27.6 --twist -166.7 --cyl_outer_diameter 200 --nr_asu 13 --z_percentage 0.3


%% Filter template of rln_104/Refine3D/job006 to 30 A resolution

% Apply lowpass filter
relion_image_handler --i /home/Medalia/Projects7/Bruno/ActinPolarity/template_library/run_class001_bundles_job017_104_mirrored_sym.mrc --o /home/Medalia/Projects7/Bruno/ActinPolarity/template_library/run_class001_bundles_job017_104_mirrored_sym_lowpass_30A.mrc --angpix 3.443 --lowpass 30

