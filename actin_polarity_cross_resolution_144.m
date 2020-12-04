

%% Actin polarity cross-resolution measurement / Box size 144 pixel

% The docking was done at choosen threshold noted below. Several click docking with
% correlation metric. The final docking ccc is also noted below.

% Chimera resampling command = vop resample #1 onGrid #0

% Docking reference = emd_6179_144_pixels_angpix_3_4A_res_8A_sym.mrc, Threshold @ 0.02

% Cross-resolution reference = emd_6179_144_pixels_angpix_3_4A_res_6_886A_sym.mrc

% Cross-resolution experiment (01) --- Resolution measurement / rln_144 / Refine3D / job006 / Template at 30 A
% Input = run_class001_bundles_job017_144_mirrored_sym_lowpass_30A.mrc, Threshold @ 0.02, ccc = 0.8361, cross_resolution_144_input_1.mrc

% Cross-resolution experiment (02) --- Resolution measurement / rln_144 / Class3D / job010
% Input = run_it050_class001.mrc, Threshold @ 0.024, ccc = 0.8844, cross_resolution_144_input_2.mrc

% Cross-resolution experiment (03) --- Resolution measurement / rln_144 / Class3D / job010 / with on purpose docked with reversed direction
% Input = run_it050_class001.mrc, Threshold @ 0.024, ccc = 0.7314, cross_resolution_144_input_3.mrc

% Cross-resolution experiment (04) --- Resolution measurement / rln_144 / Refine3D / job006
% Input = run_half1_class001_unfil.mrc, Threshold @ 0.02, ccc = 0.8701, cross_resolution_144_input_4.mrc

% Cross-resolution experiment (05) --- Resolution measurement / rln_144 / Refine3D / job006 / with on purpose docked with reversed direction
% Input = run_half1_class001_unfil.mrc, Threshold @ 0.02, ccc = 0.7292, cross_resolution_144_input_5.mrc

% Cross-resolution experiment (06) --- Resolution measurement / rln_144 / Refine3D / job013 / particles_job004_psi_prior_rev.star
% Input = run_half1_class001_unfil.mrc, Threshold @ 0.027, ccc = 0.8936, cross_resolution_144_input_6.mrc [docking with overlap metric]

% Cross-resolution experiment (07) --- Resolution measurement / rln_144 / Refine3D / job013 / particles_job004_psi_prior_rev.star
% Input = run_half1_class001_unfil.mrc, Threshold @ 0.027, ccc = 0.8862, cross_resolution_144_input_7.mrc [docking with correlation metric]

% Cross-resolution experiment (08) --- Resolution measurement / rln_144 / Refine3D / job014 / particles_job004_psi_prior_flip_rev.star
% Input = run_half1_class001_unfil.mrc, Threshold @ 0.022, ccc = 0.8678, cross_resolution_144_input_8.mrc

% Cross-resolution experiment (09) --- Resolution measurement / rln_144 / Refine3D / job015 / particles_job004_psi_prior_rev.star
% with local helical search, finer sampling and strict highres
% Input = run_half1_class001_unfil.mrc, Threshold @ 0.022, ccc = 0.883, cross_resolution_144_input_9.mrc

% Cross-resolution experiment (10) --- Resolution measurement / rln_144 / Refine3D / job017 / particles_job004_prior_particles.star
% with local helical search, finer sampling and strict highres
% Input = run_half1_class001_unfil.mrc, Threshold @ 0.022, ccc = 0.8874, cross_resolution_144_input_10.mrc

% Cross-resolution experiment (11) --- Resolution measurement / rln_144 / Refine3D / job015 / rln_postprocess / particles_job004_psi_prior_rev.star
% with local helical search, finer sampling and strict highres
% Input = postprocess_adhoc_job015_bfac.mrc, Threshold @ 0.031, ccc = 0.9024, cross_resolution_144_input_11.mrc

% Cross-resolution experiment (12) --- Resolution measurement / rln_144 / Refine3D / job017 / rln_postprocess / particles_job004_prior_particles.star
% with local helical search, finer sampling and strict highres
% Input = postprocess_adhoc_job017_bfac.mrc, Threshold @ 0.031, ccc = 0.908, cross_resolution_144_input_12.mrc

% Cross-resolution experiment (13) --- Resolution measurement / rln_144 / Refine3D / job029 / rln_postprocess / particles_job004_prior_particles_del_translational_priors_rand_red_9931.star
% with local helical search, finer sampling and strict highres (particle number restricted to same number of particles than Refine3D/job015)
% Input = postprocess_adhoc_job029_bfac.mrc, Threshold @ 0.031, ccc = 0.8985, cross_resolution_144_input_13.mrc


%% Calculate cross-resolution

% Load reference structure
ref = tom_mrcread('/home/Medalia/Projects7/Bruno/ActinPolarity/template_library/emd_6179_144_pixels_angpix_3_4A_res_6_886A_sym.mrc');
ref = single(ref.Value);

% Load (docked) result map
%result = tom_mrcread('/home/Medalia/Projects7/Bruno/ActinPolarity/cross_resolution/cross_resolution_144/cross_resolution_144_input_1.mrc');
%result = tom_mrcread('/home/Medalia/Projects7/Bruno/ActinPolarity/cross_resolution/cross_resolution_144/cross_resolution_144_input_5.mrc');
%result = tom_mrcread('/home/Medalia/Projects7/Bruno/ActinPolarity/cross_resolution/cross_resolution_144/cross_resolution_144_input_4.mrc');
%result = tom_mrcread('/home/Medalia/Projects7/Bruno/ActinPolarity/cross_resolution/cross_resolution_144/cross_resolution_144_input_6.mrc');
%result = tom_mrcread('/home/Medalia/Projects7/Bruno/ActinPolarity/cross_resolution/cross_resolution_144/cross_resolution_144_input_7.mrc');
%result = tom_mrcread('/home/Medalia/Projects7/Bruno/ActinPolarity/cross_resolution/cross_resolution_144/cross_resolution_144_input_8.mrc');
%result = tom_mrcread('/home/Medalia/Projects7/Bruno/ActinPolarity/cross_resolution/cross_resolution_144/cross_resolution_144_input_9.mrc');
%result = tom_mrcread('/home/Medalia/Projects7/Bruno/ActinPolarity/cross_resolution/cross_resolution_144/cross_resolution_144_input_10.mrc');
result = tom_mrcread('/home/Medalia/Projects7/Bruno/ActinPolarity/cross_resolution/cross_resolution_144/cross_resolution_144_input_11.mrc');%PLOT
%result = tom_mrcread('/home/Medalia/Projects7/Bruno/ActinPolarity/cross_resolution/cross_resolution_144/cross_resolution_144_input_13.mrc');%PLOT
result = single(result.Value);

% Define FSC mask / spherical
mask_fsc_sp = tom_spheremask(ones(144,144,144),52,8);

% Define FSC mask / cylindrical
mask_fsc_cyl = tom_cylindermask(ones(144,144,144),8,8);

% Normalize both maps within cylindrical mask
ref_norm = tom_norm(ref,'mean0+1std',mask_fsc_cyl);
result_norm = tom_norm(result,'mean0+1std',mask_fsc_cyl);

% Check
tom_dev(ref_norm,'',mask_fsc_cyl)
tom_dev(result_norm,'',mask_fsc_cyl)

% Create random volumes
rand_1 = (rand(144,144,144)-0.5).*10;
rand_2 = (rand(144,144,144)-0.5).*10;

% Normalize both random maps within cylindrical mask
rand_1_norm = tom_norm(rand_1,'mean0+1std',mask_fsc_cyl);
rand_2_norm = tom_norm(rand_2,'mean0+1std',mask_fsc_cyl);

% Check
tom_dev(rand_1_norm,'',mask_fsc_cyl)
tom_dev(rand_2_norm,'',mask_fsc_cyl)

% Calculate FSC
fsc = tom_compare(ref_norm.*mask_fsc_cyl.*mask_fsc_sp,result_norm.*mask_fsc_cyl.*mask_fsc_sp,72);
fsc = fsc(:,9);
fsc(1,1) = 1;

% Calculate FSC / rand control
fsc_rand = tom_compare(rand_1_norm.*mask_fsc_cyl.*mask_fsc_sp,rand_2_norm.*mask_fsc_cyl.*mask_fsc_sp,72);
fsc_rand = fsc_rand(:,9);
fsc_rand(1,1) = 1;

% Define objectpixelsize
Objectpixelsize = 3.443;%A

% Calculate resolution
f = (1:72)./(2.*Objectpixelsize.*72);
x = round((1./f).*10)./10;%A

% Add spatial frequency values
figure(1);set(gca,'XTick',18:18:72);set(gca,'XTickLabel',{'1/27.5'   '1/13.8'    '1/9.2'    '1/6.9'});

% Define 0.14 threshold criterion
t_0_14 = 0.14.*ones(1,72);

% Define 0.5 threshold criterion
t_0_5 = 0.5.*ones(1,72);

% Initialize FSC plot
figure(1);hold on;
plot(fsc,'b-','LineWidth',1);
%plot(fsc_rand,'k--');
plot(t_0_5,'k--');
plot(t_0_14,'k--');
plot(zeros(1,72),'k--');
xlim([1 72]);
ylim([-0.1 1.1]);
xlabel('Spatial frequency (1/A)');
ylabel('Correlation coefficient');
box on;

