function [omega] = stream_mod_sulq_variance(d, c, e_p, delta, scaler)
%STREAM_MOD_SULQ_VARIANCE Generate variance (β) for streaming MOD-SuLQ.
% 
% Function that generates the correct variance (omega) to be used in the 
% perturbation mask. There is also the option to set a variance scaler if
% desired.
%
% Based on work of Grammenos et al.: https://arxiv.org/abs/1907.08059
%
% Author: Andreas Grammenos (ag926@cl.cam.ac.uk)
%
% Last touched date: 19/04/2020
% 
% License: GPLv3
%

  % check if we have a specific variance scaler
  if nargin < 5
    scaler = 1;
  end

  % part of the omega variance
  part1 = (4 * d) / (e_p * c);
  part2 = sqrt( (2 * log( (d^2) / (delta * sqrt(2*pi)) )) );
  part3 = sqrt(2) / (sqrt(e_p) * c);
  
  % generate the omerga variance value
  omega = scaler * ((part1 * part2) + part3);
end

