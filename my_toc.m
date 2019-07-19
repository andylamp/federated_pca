function [t] = my_toc(ts, print, use_min)
%MY_TOC prints a formatted toc from a tic context in a similar format
% like the other messages we have. Also supports displaying only in seconds
% minutes & seconds.
%
% Author: Andreas Grammenos (ag926@cl.cam.ac.uk)
%
% Last touched date: 17/07/2019
% 
% License: GPLv3
%
  % check if the print flag is up
  if nargin < 2
    print = 1;
  end

  % setup the defaults
  if nargin < 3
    use_min = 1;
  end
  
  % grab the time from the `ts` context
  t = toc(ts);
  % check if we enabled printing
  if print == 1
    % check if we use only seconds or minutes and seconds.
    if use_min == 1
      fprintf(" -- Elapsed time %d minutes %f seconds\n", ...
        floor(t/60), rem(t, 60));
    else
      fprintf(" -- Elapsed time %d seconds..\n", t);
    end
  end
end

