function [beta_var] = mod_sulq_variance(d, n, e_p, delta)
%MOD_SULQ_VARIANCE Generate mod-sulq variance for a given setting.
%
% Based on work of Chaudhuri et al.: https://arxiv.org/pdf/1207.2812
%
% Author: Andreas Grammenos (ag926@cl.cam.ac.uk)
%
% Last touched date: 03/10/2019
% 
% License: GPLv3
%
  % variance parts as per paper
  part_a = (d + 1) / (n * e_p);
  part_b = sqrt(2 * log( (d^2 + d) / (delta * 2 * sqrt(2*pi)) ));
  part_c = 1 / (n * sqrt(e_p));
  % assign the variance as our output
  beta_var = (part_a * part_b) + part_c;
end