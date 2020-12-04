function [particle_identifier] = actin_polarity_construct_particle_identifier_prev(k,i,zaehler)
%%%%%%%%%%%%
%%%%%%%%%%%%
% This function constructs an particle identifier.
%
% INPUT
% k --- Tomogram index within BasePath
% i --- Particle index within tomogram
% zaehler --- Running particle number
%
% OUTPUT
% particle_identifier --- Name of particle

% Construct particle identifier

% Tomogram zeros
if k < 10
   tzeros = '000';
else
   if k < 100
        tzeros = '00';
   else
        if k < 1000
             tzeros = '0';
        else
             tzeros = [];
        end
   end
end

% Particle zeros
if i < 10
   pzeros = '000';
else
   if i < 100
        pzeros = '00';
   else
        if i < 1000
             pzeros = '0';
        else
             pzeros = [];
        end
   end
end

% Particle number zeros
if zaehler<10
   nzeros = '0000';
else
   if zaehler < 100
        nzeros = '000';
   else
        if zaehler < 1000
             nzeros = '00';
        else
             if zaehler < 10000
                   nzeros = '0';
             else
                   nzeros = [];
             end
        end
   end
end

particle_identifier = ['actin_particle_t_' tzeros num2str(k) '_p_' pzeros num2str(i) '_n_' nzeros num2str(zaehler)];

