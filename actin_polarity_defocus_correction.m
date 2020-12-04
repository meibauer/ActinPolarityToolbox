
% Generate BasePath
actin_polarity_construct_base_path;

% Define microscope parameter
Objectpixelsize = 0.3443;%nm
Voltage = 300;%kV
Cs = 2.7;%mm

% Define gradient convention
GradCon = +1;

% Perform CTF correction
for k=1:size(BasePath,2)
    
    % Load tiltangles and tiltaxis
    load([BasePath{k} '/alg/alg.mat'],'alg0');
    tiltangles = alg0.Tiltangles;
    tiltaxis = alg0.Tiltaxis;
    
    % Load defocus correction values / tom
    load([BasePath{k} '/ps/dz_eval_tom_mod.mat'],'dz_eval_proj_mod');
    
    % Prepare defocus parameter
    dz_1 = dz_eval_proj_mod;
    dz_2 = dz_eval_proj_mod;
    alpha = zeros(1,size(tiltangles,2));
    
    % Correct tiltseries
    CorrectTiltseries(dz_1,dz_2,alpha,tiltangles,tiltaxis,GradCon,[BasePath{k} '/proj/sorted'],'sorted_',[BasePath{k} '/proj/ctfcor'],'ctfcor_',32,512,Objectpixelsize,Voltage,Cs,4096);
    
    clear alg0;
    clear tiltangles;
    clear tiltaxis;
    clear dz_eval_proj;
    clear dz_1;
    clear dz_2;
    clear alpha;
    
    % Save ctf correction done reporter
    ctfcor = 1;
    save([BasePath{k} '/ctfcor.mat'],'ctfcor');
    clear ctfcor;
    
    disp(k);
    
end

