function [Usg, Vsg, opt_out] = my_grouse(Y, r, params)
%MY_GROUSE This is a essentially a wrapper function for GROUSE.
% It enables invoking GROUSE with the specified parameters. 
%
% GROUSE (https://arxiv.org/pdf/1702.01005.pdf) is based on work by Balzano
% et al.
%
% Author: Andreas Grammenos (ag926@cl.cam.ac.uk)
%
% Last touched date: 17/04/2020
% 
% License: GPLv3
%
    fprintf('\n ** Running GROUSE...');
    
    % default is, no error flag disabled
    if nargin < 3
      params.no_err = 1;
      params.max_rank = r;
      params.max_cycles = 1;
      params.step_size = 2;
      params.use_blk_err = 1;
    % params was supplied, check what values it has and fill in the
    % others with defauls.
    else
      % check if no err flag exists
      if ~isfield(params, 'no_err')
        params.no_err = 1;
      end
      
      % check if max rank exists
      if ~isfield(params, 'max_rank')
        params.max_rank = r;
      end
      
      % check of max cycles exists
      if ~isfield(params, 'max_cycles')
        params.max_cycles = 1;
      end
      
      % default is stepSize = 2.
      %
      % Note: in original code is 0.1 but it had worse performance.
      if ~isfield(params, 'step_size')
        params.step_size = 2;
      end
    end
    
    % check the error block size, if not present
    if ~isfield(params, 'use_blk_err')
      params.use_blk_err = 1;
    end
    
    % check the error block size, if not present
    if ~isfield(params, 'err_blk_size')
      params.err_blk_size = 100;
    end
    
    % Number of rows and columns
    [numr, numc] = size(Y);

    % grouse params
    params.en_sparse = 0;  % don't use sparse matrices
    params.rperm = 0;      % enable/disable random column permutation

    % run grouse
    ts = tic;
    [Usg, Vsg, ~, opt_out.T, opt_out.ErrFro] = grouse(1, 1, Y, numr, numc, params);
    
    % current trial execution time
    opt_out.t = my_toc(ts);
    
    fprintf("\n Finished running GROUSE ...");
end

