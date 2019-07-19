function [ T, Err, Usg, Vsg, t ] = my_grouse(Y, r, no_err, ...
  max_rank, max_cycles, step_size)
%MY_GROUSE This is a essentially a wrapper function which calls GROUSE 
% (https://arxiv.org/pdf/1702.01005.pdf) with the specified parameters 
%
% Author: Andreas Grammenos (ag926@cl.cam.ac.uk)
%
% Last touched date: 30/12/2018
% 
% License: GPLv3
%
    fprintf('\n ** Running GROUSE...');
    
    % default is, no error flag disabled
    if nargin < 3
      no_err = 0;
    end

    % default is maxRank = r
    if nargin < 4
        max_rank = r;
    end

    % default is maxCycles = 1
    if nargin < 5
        max_cycles = 1;
    end

    % default is stepSize = 2 (in original code = 0.1)
    if nargin < 6
        step_size = 2;
    end
    
    % Number of rows and columns
    [numr, numc] = size(Y);

    % grouse params
    en_sparse = 0;  % don't use sparse matrices
    rperm = 0;      % enable/disable random column permutation

    % run grouse
    ts = tic;
    [Usg, Vsg, ~, T, Err] = grouse(1, 1, Y, numr, numc, max_rank, ...
        step_size, max_cycles, en_sparse, rperm, no_err);
    t = my_toc(ts);  % current trial execution time
    fprintf("\n");
end

