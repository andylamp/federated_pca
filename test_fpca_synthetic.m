% This script is responsible for evaluating the performance of F-PCA
% against various methods when using synthetic datasets. We elected to use
% a dataset that is generated using a power law spectrum with different
% alpha values (large alpha means a very concentrated spectrum whereas a
% small alpha means a more uniform spectrum).
%
% This allows us to benchmark how our method behaves across different
% scenarios. Concretely, we compare against:
%
% - PM: https://arxiv.org/pdf/1307.0032.pdf
% - GROUSE: https://arxiv.org/pdf/1702.01005.pdf
% - SPIRIT: https://dl.acm.org/doi/10.5555/1083592.1083674
% - FD: https://arxiv.org/abs/1501.01711.pdf
%
% Based on work of Grammenos et al.: https://arxiv.org/abs/1907.08059
%
% Author: Andreas Grammenos (ag926@cl.cam.ac.uk)
%
% Last touched date: 03/06/2020
% 
% License: GPLv3
%

%% Initialisation
clc; clear; close all;

% for reproducibility
rng(300);

% use block error
use_blk_err = 0;
% enable printing
pflag = 0;
% use print figs
fig_print = 1;
% print pdfs
pdf_print = 1;

% target rank (or seed for adaptive)
rank = 5;
% number of features in each vector
feats = 400; % 1k, 200
% number of feature vectors to process
T = 4000; % 10k, 2k
% pick the spectrum range to test against
alphas = [0.001, 0.1, 0.5, 1, 2, 3];   
% the length is the range for testing
alphas_len = size(alphas, 2);
% number of simulations
nSim = 3;

% set the real dataset
params.type = "fro";
% enable printing
params.pflag = 1;
% print pdfs
params.pdf_print = 1;

% setup the environment
params = setup_vars(params);

% synthetic dataset parameters
synth_params.spectrum_type = "pl";
synth_params.lambda = 1;

% enable analytical error calculation
no_err_flag = 0;

% what to run
params.pm_run = 1;
params.sp_run = 1;
params.fd_run = 1;
params.gr_run = 1;
params.fpca_fixed_run = 0;

% what to print
params.fro_print = 1;
params.mse_print = 1;
params.subspace_err_print = 1;
params.rel_err_print = 1;
params.ext_fro_print = 1;
params.times_print = 1;


% production print, for shorter titles 
% used in the paper
params.prod_print = 1;

%% Trial variable preallocation

% preallocate fro errors
errFroFPCAFinal = NaN(nSim, alphas_len);
errFroFPCALoFinal = NaN(nSim, alphas_len);
errFroFPCAHiFinal = NaN(nSim, alphas_len);
errFroPMFinal = NaN(nSim, alphas_len); 
errFroSPFinal = NaN(nSim, alphas_len);
errFroFDFinal = NaN(nSim, alphas_len);
errFroGRFinal = NaN(nSim, alphas_len);

% preallocate mse errors
errMSEFPCAFinal = NaN(nSim, alphas_len);
errMSEFPCALoFinal = NaN(nSim, alphas_len);
errMSEFPCAHiFinal = NaN(nSim, alphas_len);
errMSEPMFinal = NaN(nSim, alphas_len); 
errMSESPFinal = NaN(nSim, alphas_len);
errMSEFDFinal = NaN(nSim, alphas_len);
errMSEGRFinal = NaN(nSim, alphas_len);

% preallocate subspace mse errors
subspaceTopRFPCAFinal = NaN(nSim, alphas_len);
subspaceTopRFPCALoFinal = NaN(nSim, alphas_len);
subspaceTopRFPCAHiFinal = NaN(nSim, alphas_len);
subspaceTopRPMFinal = NaN(nSim, alphas_len);  
subspaceTopRSPFinal = NaN(nSim, alphas_len);
subspaceTopRFDFinal = NaN(nSim, alphas_len);
subspaceTopRGRFinal = NaN(nSim, alphas_len);

% relative erros for SPIRIT and F-PCA over T
relerrs_fpca = NaN(alphas_len, nSim, T);
relerrs_sp = NaN(alphas_len, nSim, T);

% Fro erros for SPIRIT, PM, and F-PCA over T
fro_err_fpca = NaN(alphas_len, nSim, T);
fro_err_sp = NaN(alphas_len, nSim, T);
% fro_err_pm_lo = NaN(alphas_len, nSim, T);
% fro_fd_lo = NaN(alphas_len, nSim, T);

% Principal Component evolution counters over T
rpcs_fpca = NaN(alphas_len, nSim, T);
rpcs_sp = NaN(alphas_len, nSim, T);

% timers
t_fpca = NaN(nSim, alphas_len);
t_fpca_lo = NaN(nSim, alphas_len);
t_fpca_hi = NaN(nSim, alphas_len);
t_spirit = NaN(nSim, alphas_len);
t_pm = NaN(nSim, alphas_len);
t_fd = NaN(nSim, alphas_len);
t_gr = NaN(nSim, alphas_len);

% xt axis for figures
xt = 1:alphas_len;

% figure legends pre-process based on what we run
legendCells = {'SP', 'PM', 'FD', 'GROUSE', ...
  'F-PCA', 'F-PCA_{lo}', 'F-PCA_{hi}'};

% check if we run the power method
if params.pm_run == 0
  idc = ismember(legendCells, {'PM'});
  legendCells = legendCells(~idc);
end

% check if we run frequent directions
if params.fd_run == 0
  idc = ismember(legendCells, {'FD'});
  legendCells = legendCells(~idc);
end

% check if we run 
if params.gr_run == 0
  idc = ismember(legendCells, {'FD'});
  legendCells = legendCells(~idc);
end

if params.sp_run == 0
  idc = ismember(legendCells, {'SP'});
  legendCells = legendCells(~idc);
end


if params.fpca_fixed_run == 0
  % remove low
  idc = ismember(legendCells, {'F-PCA_{lo}'});
  legendCells = legendCells(~idc);
  % and  hi
  idc = ismember(legendCells, {'F-PCA_{hi}'});
  legendCells = legendCells(~idc);
end


%% Trial execution

for trial = 1:size(alphas, 2)
  
    fprintf("\n -- Running for alpha %d\n\n", alphas(trial));
    
  for csim = 1:nSim
    
    % Generate matrix based on spectrum
    fprintf("\n !!! Running simulation %d out of %d using alpha %d !!!\n", ...
      csim, nSim, alphas(trial));
    % assign the alpha for the current run
    synth_params.alpha = alphas(trial);
    % generate the data
    [Y, Sigma] = synthetic_data_gen(feats, T, synth_params);

    % Test F-PCA edge
    fprintf("\n\n -- Running F-PCA Edge test\n\n");

    fprintf("\n ** Running with rank seed of %d (out of: %d)\n", ...
      rank, feats);

    % F-PCA parameters
    clear fpca_params
    fpca_params.adaptive = 1;
    fpca_params.holdoff = 1;
    fpca_params.blk_size = 50;
    fpca_params.no_err = no_err_flag;
    fpca_params.use_blk_err = use_blk_err;
    fpca_params.use_ext_errs = 1;
    fpca_params.verbose = 0;

    [Ufpca, ~, fpca_opt_out] = fpca_edge(Y, rank, fpca_params);

    fro_err_fpca(trial, csim, :) = fpca_opt_out.ErrFro;
    t_fpca(csim, trial) = fpca_opt_out.t;
    rmax = fpca_opt_out.rmax;

    relerrs_fpca(trial, csim, :) = fpca_opt_out.relerrs;
    rpcs_fpca(trial, csim, :) = fpca_opt_out.rpcs;

    % compute the fro and mse errors
    yr_sz = size(fpca_opt_out.Yr, 2);
    yr_err_sum = sum(sum((Y(:, 1:yr_sz)-fpca_opt_out.Yr).^2, 1));
    fpca_fro = yr_err_sum / yr_sz;
    fpca_mse = yr_sz * immse(Y(:, 1:yr_sz), fpca_opt_out.Yr);
    fprintf(" !! Final fro: %d (with seed %d)\n", ...
      fpca_fro, rank);

    % tag the errors
    errFroFPCAFinal(csim, trial) = fpca_fro;
    errMSEFPCAFinal(csim, trial) =  fpca_mse;

    fprintf("\n\n -- Finished running F-PCA edge test\n");

    % Test F-PCA
    if params.fpca_fixed_run == 1

      fprintf(" -- Running Fixed F-PCA Edge test\n\n");

      lo_r = min(rank, size(Ufpca, 2));
      fprintf("\n ** Running low rank (r: %d)\n", lo_r); 
      
      clear fpca_lo_params
      fpca_lo_params.adaptive = 0;
      fpca_lo_params.holdoff = 1;
      fpca_lo_params.blk_size = 50;
      fpca_lo_params.no_err = no_err_flag;
      fpca_lo_params.use_ext_errs = 1;
      fpca_lo_params.use_blk_err = use_blk_err;
      fpca_lo_params.verbose = 0;

      [Ufpca_lo, ~, fpca_lo_opt_out] = fpca_edge(Y, lo_r, fpca_lo_params);

      lo_yr = size(fpca_lo_opt_out.Yr, 2);
      fro_fpca_lo = sum(sum((Y(:, 1:lo_yr) - fpca_lo_opt_out.Yr).^2, 1)) / lo_yr;
      mse_fpca_lo = lo_yr * immse(Y(:, 1:lo_yr), fpca_lo_opt_out.Yr);
      fprintf(" !! Final fro: %d (low_r: %d)\n", fro_fpca_lo, lo_r);
      
      % get the time
      t_fpca_lo(csim, trial) = fpca_lo_opt_out.t;

      errFroFPCALoFinal(csim, trial) = fro_fpca_lo;
      errMSEFPCALoFinal(csim, trial) = mse_fpca_lo;

      % need to run high rank to bound it, so we 
      % set the high rank to be the max rank observed
      hi_r = rmax;

      fprintf("\n ** Running high rank (r: %d) to bound\n", rmax);

      clear fpca_hi_params
      fpca_hi_params.adaptive = 0;
      fpca_hi_params.holdoff = 1;
      fpca_hi_params.blk_size = 50;
      fpca_hi_params.no_err = no_err_flag;
      fpca_hi_params.use_blk_err = use_blk_err;
      fpca_hi_params.use_ext_errs = 1;
      fpca_hi_params.verbose = 0;

      [Ufpca_hi, ~, fpca_hi_opt_out] = fpca_edge(Y, hi_r, fpca_hi_params);

      % get the size
      hi_yr = size(fpca_hi_opt_out.Yr, 2);

      % get the time
      t_fpca_hi(csim, trial) = fpca_hi_opt_out.t;

      % compute the frobenious error
      fro_fpca_hi = sum(sum((Y(:, 1:hi_yr) - fpca_hi_opt_out.Yr).^2, 1)) / hi_yr;
      fprintf(" !! Final fro: %d (rank: %d)\n", fro_fpca_hi, hi_r);

      % compute the mse
      mse_fpca_hi = hi_yr * immse(Y(:, 1:hi_yr), fpca_hi_opt_out.Yr);

      errFroFPCAHiFinal(csim, trial) = fro_fpca_hi;
      errMSEFPCAHiFinal(csim, trial) = mse_fpca_hi;


      fprintf("\n -- Finished running Fixed F-PCA edge test\n");
    end

    % Test Mitliagkas Power Method
    if params.pm_run == 1
      % pm_rank = % min(pm_rank, lo_r);
      fprintf("\n -- Running PM with rank %d (out of: %d)\n", ...
        rank, feats);

      clear pm_params;
      pm_params.no_err = no_err_flag;
      pm_params.use_blk_err = use_blk_err;
      
      [Upms, pm_opt_out] = mitliag_pm(Y, rank, pm_params);
      
      % get the time
      t_pm(csim, trial) = pm_opt_out.t;

      lo_yr = size(pm_opt_out.Yr, 2);
      fro_pm = sum(sum((Y(:, 1:lo_yr)-pm_opt_out.Yr).^2, 1)) / lo_yr;
      mse_pm = lo_yr * immse(Y(:, 1:lo_yr), pm_opt_out.Yr);

      errFroPMFinal(csim, trial) = fro_pm;
      errMSEPMFinal(csim, trial) = mse_pm;
      fprintf(" !! Final fro: %d (rank: %d)\n", fro_pm, rank);
      fprintf("\n -- Finished Running PM\n");
    end
    
    % Test FD
    if params.fd_run == 1
      fprintf("\n -- Running FD with rank %d (out of %d)\n", ...
        rank, feats);
      % run fd
      clear fd_params;
      fd_params.no_err = 0;
      fd_params.use_blk_err = use_blk_err;
      
      [Ufd, fd_opt_out] = fd(Y', rank, fd_params);
      % transpose it
      fd_opt_out.Yr = fd_opt_out.Yr';
      Ufd = Ufd';
      
      % get the fd time
      t_fd(csim, trial) = fd_opt_out.t;
      
      fd_yr = size(fd_opt_out.Yr, 2);
      fro_fd = sum(sum((Y(:, 1:fd_yr)-fd_opt_out.Yr).^2, 1))/fd_yr;
      mse_fd = fd_yr * immse(Y(:, 1:fd_yr), fd_opt_out.Yr);
      
      errFroFDFinal(csim, trial) = fro_fd;
      errMSEFDFinal(csim, trial) = mse_fd;
      fprintf(" !! Final fro: %d (rank: %d)\n", fro_fd, rank);
      fprintf("\n -- Finished Running FD\n");
    end
    
    % Test GROUSE
    if params.gr_run == 1
      fprintf("\n -- Running GROUSE with rank %d (out of %d)\n", ...
        rank, feats);
      
      % run grouse
      clear gr_params;
      gr_params.no_err = no_err_flag;
      gr_params.use_blk_err = use_blk_err;
      
      [U_gr, V_gr, gr_opt_out] = my_grouse(Y, rank, gr_params);
      
      t_gr(csim, trial) = gr_opt_out.t;
      
      % expand U_gr*V_gr' to get the Yr_gr
      Yr_gr = U_gr*V_gr';   
      gr_yr = size(Yr_gr, 2);
      % now get the fro error
      fro_gr_lo_r = sum(sum((Y(:, 1:gr_yr)-Yr_gr).^2, 1)) / gr_yr;
      mse_gr = gr_yr * immse(Y(:, 1:gr_yr), Yr_gr); 
      
      errFroGRFinal(csim, trial) = fro_gr_lo_r;
      errMSEGRFinal(csim, trial) = mse_gr;
      
      fprintf(" !! Final fro: %d (rank %d)\n", fro_gr_lo_r, rank);
      fprintf("\n -- Finished Running GROUSE\n");
    end


    % Test SPIRIT
    if params.sp_run == 1
      fprintf("\n -- Running SPIRIT test\n\n");
      
      % SPIRIT parameters
      % sp_lambda = .9; % use 1.0 to exhibit pathological behaviour
      sp_lambda = 1.0;
      sp_energy = [0.95, 0.98];
      
      clear sp_params;
      sp_params.no_err = no_err_flag;
      sp_params.use_blk_err = use_blk_err;
      
      [W1, sp_opt_out] = SPIRIT(Y', sp_lambda, sp_energy, sp_params);
      
      % get the final run
      sp_r = sp_opt_out.final_rank;
      
      % assign fro errors
      fro_err_sp(trial, csim, :) = sp_opt_out.ErrFro;
      % assign the related errors
      relerrs_sp(trial, csim, :) = sp_opt_out.relErrors;
      % assign rpcs
      rpcs_sp(trial, csim, :) = sp_opt_out.pcs;
      
      
      % assign the ticking
      t_spirit(csim, trial) = sp_opt_out.t; 
      
      % error calculation
      YSpiritSubRecon = (W1*W1')*Y;
      fro_sp = sum(sum((Y-YSpiritSubRecon).^2, 1))/T;
      errFroSPFinal(csim, trial) = fro_sp;
      errMSESPFinal(csim, trial) = T * immse(Y, YSpiritSubRecon);
      fprintf(" !! Final fro: %d (final rank: %d)\n", fro_sp, sp_r);
      fprintf("\n -- Finished Running SPIRIT\n");
    end

    %% Test the subspaces MSE

    if params.subspace_err_print == 1

      % Compute the offline Subspace of Y
      [Uoff, ~, ~] = svd(Y);
      % expand it, for comparison against the other subspaces
      pcaUU = Uoff*Uoff';
      % to have the absolute value
      Uoff_abs = abs(Uoff);
      % Subspace for SPIRIT
      if params.sp_run == 1
        %subTopRMSE = mse(pcaUU, W1*W1');
        subTopRMSE = immse(Uoff_abs(:, 1:sp_r), abs(W1));
        fprintf("\n ** SPIRIT Subspace (for r: %d) MSE: %d", ...
          sp_r, subTopRMSE);
        subspaceTopRSPFinal(csim, trial) = subTopRMSE;
      end

      % Subspace for Power Method
      if params.pm_run == 1
        [~, r] = size(Upms);
        subTopRMSE = immse(Uoff_abs(:, 1:r), abs(Upms));
        fprintf("\n ** Mitliagkas Subspace (r: %d) MSE: %d", ...
          r, subTopRMSE);
        subspaceTopRPMFinal(csim, trial) = subTopRMSE;
      end
      
      % Subspace for FD
      if params.fd_run == 1
        [~, r] = size(Ufd);
        %subTopRMSE = mse(pcaUU, Ufd*Ufd');
        subTopRMSE = immse(Uoff_abs(:, 1:r), abs(Ufd));
        fprintf("\n ** FD Subspace (r: %d) MSE: %d", ...
          r, subTopRMSE);
        subspaceTopRFDFinal(csim, trial) = subTopRMSE;
      end
      
      % Subspace for GROUSE
      if params.gr_run == 1
        [~, r] = size(U_gr);
        subTopRMSE = immse(Uoff_abs(:, 1:r), abs(U_gr));
        fprintf("\n ** GROUSE Subspace (r: %d) MSE: %d", ...
          r, subTopRMSE);
        subspaceTopRGRFinal(csim, trial) = subTopRMSE;
      end

      % Subspace for SPCA Low
      if params.fpca_fixed_run == 1
        % lower bound
        [~, r] = size(Ufpca_lo);
        subTopRMSE = immse(Uoff_abs(:, 1:r), abs(Ufpca_lo));
        %subTopRMSE = mse(pcaUU, Ufm_lo*Ufm_lo');
        fprintf("\n ** F-PCA Subspace (lo_r: %d) MSE: %d", ...
          r, subTopRMSE);
        subspaceTopRFPCALoFinal(csim, trial) = subTopRMSE;
        % higher bound
        [~, r] = size(Ufpca_hi);
        subTopRMSE = immse(Uoff_abs(:, 1:r), abs(Ufpca_hi));
        fprintf("\n ** F-PCA Subspace UU' (hi_r: %d) MSE: %d", ...
          r, subTopRMSE);
        subspaceTopRFPCAHiFinal(csim, trial) = subTopRMSE;

      end

      % Subspace for F-PCA
      [~, r] = size(Ufpca);
      subTopRMSE = immse(Uoff_abs(:, 1:r), abs(Ufpca));
      fprintf("\n ** F-PCA Subspace (r: %d) MSE: %d", ...
        r, subTopRMSE);  
      subspaceTopRFPCAFinal(csim, trial) = subTopRMSE;

    end
    % to end with a nice console offset
    fprintf("\n");
  
  end

end

%% Plot frobenious norm errors


% print the Yr vs Y fro final error
if params.fro_print == 1
  logErrFroSPFinal = log(errFroSPFinal);
  logErrFroFPCAFinal = log(errFroFPCAFinal);
  logErrFroFPCALoFinal = log(errFroFPCALoFinal);
  logErrFroFPCAHiFinal = log(errFroFPCAHiFinal);
  logErrFroGRFinal = log(errFroGRFinal);
  logErrFroPMFinal = log(errFroPMFinal);
  logErrFroFDFinal = log(errFroFDFinal);

  mlogErrFroSPFinal = mean(logErrFroSPFinal, 1);
  mlogErrFroFPCAFinal = mean(logErrFroFPCAFinal, 1);
  mlogErrFroFPCALoFinal = mean(logErrFroFPCALoFinal, 1);
  mlogErrFroFPCAHiFinal = mean(logErrFroFPCAHiFinal, 1);
  mlogErrFroGRFinal = mean(logErrFroGRFinal, 1);
  mlogErrFroPMFinal = mean(logErrFroPMFinal, 1);
  mlogErrFroFDFinal = mean(logErrFroFDFinal, 1);

  slogErrFroSPFinal = std(logErrFroSPFinal, 1);
  slogErrFroFPCAFinal = std(logErrFroFPCAFinal, 1);
  slogErrFroFPCALoFinal = std(logErrFroFPCALoFinal, 1);
  slogErrFroFPCAHiFinal = std(logErrFroFPCAHiFinal, 1);
  slogErrFroGRFinal = std(logErrFroGRFinal, 1);
  slogErrFroPMFinal = std(logErrFroPMFinal, 1);
  slogErrFroFDFinal = std(logErrFroFDFinal, 1);
  
  fig = figure;
  hold on;
  
  if params.sp_run == 1
    errorbar(xt, mlogErrFroSPFinal, slogErrFroSPFinal, '--');
  end
  
  if params.pm_run == 1
    errorbar(xt, mlogErrFroPMFinal, slogErrFroPMFinal);
  end
  
  if params.fd_run == 1
    errorbar(xt, mlogErrFroFDFinal, slogErrFroFDFinal);
  end
  
  if params.gr_run == 1
    errorbar(xt, mlogErrFroGRFinal, slogErrFroGRFinal);
  end
  
  errorbar(xt, mlogErrFroFPCAFinal, slogErrFroFPCAFinal, 'LineWidth', 2);
  
  if params.fpca_fixed_run == 1
    errorbar(xt, mlogErrFroFPCALoFinal, slogErrFroFPCALoFinal);
    errorbar(xt, mlogErrFroFPCAHiFinal, slogErrFroFPCAHiFinal);
  end
  hold off;
  
  if params.prod_print == 0
    stl = sprintf(['alphas vs fro errors ' , ...
      '(seed rank: %d, sims: %d, T: %d, feats: %d)'], ...
      rank, nSim, T, feats);
  else
    stl = sprintf('fro errors for varying alpha');
  end
  title(stl);
  xticks(1:alphas_len);
  xticklabels(alphas);
  ylabel('fro error');
  xlabel('alphas');
  
  % put legends
  legend(legendCells, 'Location', 'best');
  
  st = sprintf("fro_err_rseed_%d_sims_%d_T_%d", rank, nSim, T);
  print_fig(fig, st, params);
end

%% Plot MSE errors

if params.mse_print == 1
  % take the log scale of the measurements
  logErrMSESPFinal = log(errMSESPFinal);
  logErrMSEFPCAFinal = log(errMSEFPCAFinal);
  logErrMSEFPCALoFinal = log(errMSEFPCALoFinal);
  logErrMSEFPCAHiFinal = log(errMSEFPCAHiFinal);
  logErrMSEPMFinal = log(errMSEPMFinal);
  logErrMSEFDFinal = log(errMSEFDFinal);
  logErrMSEGRFinal = log(errMSEGRFinal);

  % generate the mean
  mlogErrMSESPFinal = mean(logErrMSESPFinal, 1);
  mlogErrMSEFPCAFinal = mean(logErrMSEFPCAFinal, 1);
  mlogErrMSEFPCALoFinal = mean(logErrMSEFPCALoFinal, 1);
  mlogErrMSEFPCAHiFinal = mean(logErrMSEFPCAHiFinal, 1);
  mlogErrMSEPMFinal = mean(logErrMSEPMFinal, 1);
  mlogErrMSEFDFinal = mean(logErrMSEFDFinal, 1);
  mlogErrMSEGRFinal = mean(logErrMSEGRFinal, 1);
  
  % and the standard deviation
  slogErrMSESPFinal = std(logErrMSESPFinal, 1);
  slogErrMSEFPCAFinal = std(logErrMSEFPCAFinal, 1);
  slogErrMSEFPCALoFinal = std(logErrMSEFPCALoFinal, 1);
  slogErrMSEFPCAHiFinal = std(logErrMSEFPCAHiFinal, 1);
  slogErrMSEPMFinal = std(logErrMSEPMFinal, 1);
  slogErrMSEFDFinal = std(logErrMSEFDFinal, 1);
  slogErrMSEGRFinal = std(logErrMSEGRFinal, 1);

  fig = figure;
  hold on;
  
  if params.sp_run == 1
    errorbar(xt, mlogErrMSESPFinal, slogErrMSESPFinal, '--');
  end
  
  if params.pm_run == 1
    errorbar(xt, mlogErrMSEPMFinal, slogErrMSEPMFinal);
  end
  
  if params.fd_run  == 1
    errorbar(xt, mlogErrMSEFDFinal, slogErrMSEFDFinal);
  end
  
  if params.gr_run == 1
    errorbar(xt, mlogErrMSEGRFinal, slogErrMSEGRFinal);
  end
  
  errorbar(xt, mlogErrMSEFPCAFinal, slogErrMSEFPCAFinal, 'LineWidth', 2);
  
  if params.fpca_fixed_run == 1
    errorbar(xt, mlogErrMSEFPCALoFinal, slogErrMSEFPCALoFinal);
    errorbar(xt, mlogErrMSEFPCAHiFinal, slogErrMSEFPCAHiFinal);
  end
  
  hold off;
  stl = sprintf(['Yr vs real MSE ', ...
    '(seed rank: %d, sims: %d, T: %d, feats: %d)'], ...
    rank, nSim, T, feats);
  title(stl);
  xticks(1:alphas_len);
  xticklabels(alphas);
  ylabel('error (log(mse))');
  xlabel('\alpha');
  
  % put legends
  legend(legendCells, 'Location', 'best');
  
  st = sprintf("mse_err_rseed_%d_sims_%d_T_%d", rank, nSim, T);
  print_fig(fig, st, params);

end

%% Plot the subspace MSE errors

if params.subspace_err_print == 1
  % take the log scale of the measurements
  logErrsSubSPFinal = log(subspaceTopRSPFinal);
  logErrSubFPCAFinal = log(subspaceTopRFPCAFinal);
  logErrSubFPCALoFinal = log(subspaceTopRFPCALoFinal);
  logErrSubFPCAHiFinal = log(subspaceTopRFPCAHiFinal);
  logErrSubPMFinal = log(subspaceTopRPMFinal);
  logErrSubFDFinal = log(subspaceTopRFDFinal);
  logErrSubGRFinal = log(subspaceTopRGRFinal);
  
  % generate the mean
  mlogErrsSubSPFinal = mean(logErrsSubSPFinal, 1);
  mlogErrSubFPCAFinal = mean(logErrSubFPCAFinal, 1);
  mlogErrSubFPCALoFinal = mean(logErrSubFPCALoFinal, 1);
  mlogErrSubFPCAHiFinal = mean(logErrSubFPCAHiFinal, 1);
  mlogErrSubPMFinal = mean(logErrSubPMFinal, 1);
  mlogErrSubFDFinal = mean(logErrSubFDFinal, 1);
  mlogErrSubGRFinal = mean(logErrSubGRFinal, 1);
  
  % and the standard deviation
  slogErrsSubSPFinal = std(logErrsSubSPFinal, 1);
  slogErrSubFPCAFinal = std(logErrSubFPCAFinal, 1);
  slogErrSubFPCALoFinal = std(logErrSubFPCALoFinal, 1);
  slogErrSubFPCAHiFinal = std(logErrSubFPCAHiFinal, 1);
  slogErrSubPMFinal = std(logErrSubPMFinal, 1);
  slogErrSubFDFinal = std(logErrSubFDFinal, 1);
  slogErrSubGRFinal = std(logErrSubGRFinal, 1); 
  
  fig = figure;
  hold on;
  
  if params.sp_run == 1
    errorbar(xt, mlogErrsSubSPFinal, slogErrsSubSPFinal, '--');
  end
  
  if params.pm_run == 1
    errorbar(xt, mlogErrSubPMFinal, slogErrSubPMFinal);
  end
  
  if params.fd_run == 1
    errorbar(xt, mlogErrSubFDFinal, slogErrSubFDFinal);
  end
  
  if params.gr_run == 1
    errorbar(xt, mlogErrSubGRFinal, slogErrSubGRFinal);
  end
  
  errorbar(xt, mlogErrSubFPCAFinal, slogErrSubFPCAFinal, 'LineWidth', 2);
  
  if params.fpca_fixed_run == 1
    errorbar(xt, mlogErrSubFPCALoFinal, slogErrSubFPCALoFinal);
    errorbar(xt, mlogErrSubFPCAHiFinal, slogErrSubFPCAHiFinal);
  end
  
  hold off;
  if params.prod_print == 0
    stl = sprintf(['semilogy subspace vs real MSE ', ... 
      '(seed rank: %d, sims: %d, T: %d, feats: %d)'], ... 
      rank, nSim, T, feats);
  else
    stl = sprintf('Estimated vs Real Subspace MSE');
  end
  title(stl);
  xticks(1:alphas_len);
  xticklabels(alphas);
  ylabel('error (log(rmse))');
  xlabel('\alpha'); 
  
  % put legends
  legend(legendCells, 'Location', 'best');
  % print to a file
  st = sprintf("subspace_semilogy_err_feat_%d_T_%d_nsims_%d", ...
    feats, T, nSim);
  print_fig(fig, st, params);
end

%% Plot the rel. proj. errors and PC's kept over T (F-PCA vs SPIRIT only)

% rel errors print
if params.rel_err_print == 1 && params.sp_run == 1
  for i = 1:alphas_len
    fig = figure; 
    % pc evolution
    subplot(2, 1, 1);
      hold on;
      plot(mean(squeeze(rpcs_sp(i, :, :)), 1), '--', 'LineWidth', 2); 
      plot(mean(squeeze(rpcs_fpca(i, :, :)), 1), ...
         'LineWidth', 2, 'color', [0.4660    0.6740    0.1880]); 
      legend('SP', 'F-PCA', 'Location', 'best'); 
      % print to a file
      stl = sprintf(['PCs evolution over T (%d) ', ...
        'with %d feats for alpha: %6.5f'], ...
        T, feats, alphas(i));
      title(stl);
      ylabel('PC count');
      xlabel('time ticks');
      hold off;

    % rel error evolution
    subplot(2, 1, 2);
      hold on;
      plot(mean(squeeze(relerrs_sp(i, :, :)), 1), '--', 'LineWidth', 2);
      plot(mean(squeeze(relerrs_fpca(i, :, :)), 1), ...
         'LineWidth', 2, 'color', [0.4660    0.6740    0.1880]);
      legend('SP', 'F-PCA', 'Location', 'best'); 
      stl = sprintf('Rel errors over T (%d) for alpha: %6.5f', T, alphas(i));
      title(stl);
      ylabel('rel. error');
      xlabel('time ticks');
      hold off;

      % print to a file
      st = sprintf("rel_errors_feat_%d_T_%d_alpha_%s_nsims_%d", ...
        feats, T, strrep(num2str(alphas(i)), ".", "_"), nSim);
    print_fig(fig, st, params);
  end
end

%% Print rel. fro error and PC's kept over T (F-PCA vs SPIRIT only)

% pc number over time print
if params.ext_fro_print == 1 && params.sp_run == 1
  for i = 1:alphas_len
    fig = figure; 
    % pc evolution
    subplot(2, 1, 1);
      hold on;
       plot(mean(squeeze(rpcs_sp(i, :, :)), 1), '--', 'LineWidth', 2);
       plot(mean(squeeze(rpcs_fpca(i, :, :)), 1), ...
         'LineWidth', 2, 'color', [0.4660    0.6740    0.1880]); 
      
      legend('SP', 'F-PCA', 'Location', 'best'); 
      % print to a file
      stl = sprintf(['PCs evolution over T (%d) ', ...
        'with %d feats for alpha: %6.5f'], ...
        T, feats, alphas(i));
      title(stl);
      ylabel('PC count');
      xlabel('time ticks');
      hold off;

    % rel error evolution
    subplot(2, 1, 2);
      hold on;
      plot(mean(squeeze(fro_err_sp(i, :, :)), 1), '--', 'LineWidth', 2);
      plot(mean(squeeze(fro_err_fpca(i, :, :)), 1), ...
        'LineWidth', 2, 'color', [0.4660    0.6740    0.1880]);
      legend('SP', 'F-PCA', 'Location', 'best');
      stl = sprintf('Fro errs over T (%d) for alpha: %6.5f', T, alphas(i));
      title(stl);
      ylabel('error (fro)');
      xlabel('time ticks');
      hold off;

      % print to a file
      st = sprintf("fro_ext_errors_feat_%d_T_%d_alpha_%s_nsims_%d", ...
        feats, T, strrep(num2str(alphas(i)), ".", "_"), nSim);
    print_fig(fig, st, params);
  end
end

%% Print the average execution times

if params.times_print == 1
  fig = figure;
  hold on;
  
  if params.sp_run == 1
    errorbar(xt, mean(t_spirit, 1), std(t_spirit, 1));
  end
  
  if params.pm_run == 1
    errorbar(xt, mean(t_pm, 1), std(t_pm, 1));
  end
  
  if params.fd_run == 1
    errorbar(xt, mean(t_fd, 1), std(t_fd, 1));
  end
  
  if params.gr_run == 1
    errorbar(xt, mean(t_gr, 1), std(t_gr, 1));
  end
  
  errorbar(xt, mean(t_fpca, 1), std(t_fpca, 1), 'LineWidth', 2);
  
  % check if we run fpca fixed
  if params.fpca_fixed_run == 1
    errorbar(xt, mean(t_fpca_lo, 1), std(t_fpca_lo, 1));
    errorbar(xt, mean(t_fpca_hi, 1), std(t_fpca_hi, 1));
  end
  
  hold off;
  stl = sprintf("Mean Execution times (sims: %d, T: %d, feats: %d)", ...
    nSim, T, feats);
  title(stl);
  xticks(xt);
  xticklabels(alphas);
  ylabel('execution time (s)');
  xlabel('alphas');
  
  % put the legends
  legend(legendCells);
  
  % print to a file
  st = sprintf("exec_times_feat_%d_T_%d_nsims_%d", ...
    feats, T, nSim);
  print_fig(fig, st, params);
end


