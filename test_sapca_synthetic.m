%% Federated Streaming Adaptive PCA - Synthetic Dataset Tests
%
% Description: 
%   This code is supplied as additional material alongside our paper:
%   "Federated PCA with Adaptive Rank Estimation"
%   
% This script constructs a virtual federated hierarchy where the PCA 
% object is searched. The result is propageted from bottom to top and the
% final result is stored over at the root node. The result can then be used 
% to globally or selectively enhance modes located at other places by using
% our weighted subspace merging algorithm.
%
% Notes:
%  i) Please ensure you have an up-to-date MATLAB version (> 2017a) as 
%     older versions have a problems handling character/string arrays in  
%     certain cases which are extensively used in this script.
% ii) This script uses parallel pool and results may differ depending on
%     CPU & RAM resources available.
%
% Author: Andreas Grammenos (ag926@cl.cam.ac.uk)
%
% Last touched date: 17/07/2019
% 
% License: 
%  code: GPLv3, author: A. Grammenos 
%  paper: A. Grammenos, R. Mendoza-Smith, C. Mascolo, and J. Crowcroft. 
%         Authors retain their respective copyrights. 
%         
%         Pre-print link: https://arxiv.org/abs/1907.08059
%

%% Initialisation
clc; clear; close all;

% scope in globals
global use_blk_err
global pflag
global fig_print
global pdf_print

% for reproducibility
rng(300);

% use block error
use_blk_err = 0;
% enable printing
pflag = 1;
% use print figs
fig_print = 1;
% print pdfs
pdf_print = 1;

% target rank (or seed for adaptive)
r_seed = 10;
% number of features in each vector
feats = 400; % 1k, 200
% number of feature vectors to process
T = 4000; % 10k, 2k
% pick the spectrum range to test against
% alphas = [1, 2];
% alphas = [0.01, 0.5, 1, 2, 3, 4];
alphas = [0.0001, 0.001, 0.5, 1, 2, 3];   
% the length is the range for testing
alphas_len = size(alphas, 2);
% number of simulations
nSim = 5;

% enable analytical error calculation
no_err_flag = 0;

% SPIRIT parameters
sp_holdOffTime = 0;
sp_k0 = r_seed;
sp_lambda = 1; %1
%sp_lambda = 0.9; %1
sp_energy = [0.95, 0.98];
sp_silent = 1;
sp_no_err = no_err_flag;

% SAPCA parameters
am_r_seed = r_seed;
am_tr_bounds = [1, 10]; % percentages
am_holdoff = 0.9;
am_blk = 50;
am_floor_mul = 2;
am_no_err = no_err_flag;
am_silent = 1;

% SPCA Parameters
fm_rank = r_seed;
fm_blk = 50;
fm_floor_mul = 2;
fm_no_err = no_err_flag;

% GROUSE Parameters
gr_rank = r_seed;
gr_no_err = no_err_flag;

% PM Parameters
pm_rank = r_seed;
pm_blk = feats;
pm_floor_mul = 2;
pm_no_err = no_err_flag;

% FD Parameters
fd_rank = r_seed;
fd_no_err = no_err_flag;

% what to run
sapca_run = 1;
spca_run = 1;
pm_run = 1;
sp_run = 1;
fd_run = 1;
gr_run = 1;
fro_print = 1;
ext_fro_print = 1;
mse_print = 0;
subspace_err_print = 1;
rel_err_print = 1;
times_print = 1;

% production print, for shorter titles 
% used in the paper
prod_print = 1;

% setup the environment
setup_vars;

%% Trial variable preallocation

% preallocate fro errors
errFroAMFinal = NaN(nSim, alphas_len);
errFroFMLoFinal = NaN(nSim, alphas_len);
errFroFMHiFinal = NaN(nSim, alphas_len);
errFroPMFinal = NaN(nSim, alphas_len); 
errFroSPFinal = NaN(nSim, alphas_len);
errFroFDFinal = NaN(nSim, alphas_len);
errFroGRFinal = NaN(nSim, alphas_len);

% preallocate mse errors
errMSEAMFinal = NaN(nSim, alphas_len);
errMSEFMLoFinal = NaN(nSim, alphas_len);
errMSEFMHiFinal = NaN(nSim, alphas_len);
errMSEPMFinal = NaN(nSim, alphas_len); 
errMSESPFinal = NaN(nSim, alphas_len);
errMSEFDFinal = NaN(nSim, alphas_len);
errMSEGRFinal = NaN(nSim, alphas_len);

% preallocate subspace mse errors
subspaceTopRAMFinal = NaN(nSim, alphas_len);
subspaceTopRFMLoFinal = NaN(nSim, alphas_len);
subspaceTopRFMHiFinal = NaN(nSim, alphas_len);
subspaceTopRPMFinal = NaN(nSim, alphas_len);  
subspaceTopRSPFinal = NaN(nSim, alphas_len);
subspaceTopRFDFinal = NaN(nSim, alphas_len);
subspaceTopRGRFinal = NaN(nSim, alphas_len);

% relative erros for SPIRIT and SAPCA over T
relerrs_am = NaN(alphas_len, nSim, T);
relerrs_sp = NaN(alphas_len, nSim, T);

% Fro erros for SPIRIT, PM, SPCA, and SAPCA over T
fro_err_am = NaN(alphas_len, nSim, T);
fro_err_sp = NaN(alphas_len, nSim, T);
fro_err_fm_lo = NaN(alphas_len, nSim, T);
fro_err_fm_hi = NaN(alphas_len, nSim, T);
fro_err_pm_lo = NaN(alphas_len, nSim, T);
fro_fd_lo = NaN(alphas_len, nSim, T);

% Principal Component evolution counters over T
rpcs_am = NaN(alphas_len, nSim, T);
rpcs_sp = NaN(alphas_len, nSim, T);

% timers
t_sapca = NaN(nSim, alphas_len);
t_spca_lo = NaN(nSim, alphas_len);
t_spca_hi = NaN(nSim, alphas_len);
t_spirit = NaN(nSim, alphas_len);
t_pm = NaN(nSim, alphas_len);
t_fd = NaN(nSim, alphas_len);
t_gr = NaN(nSim, alphas_len);

%% Trial execution

for trial = 1:size(alphas, 2)
  
    fprintf("\n -- Running for alpha %d\n\n", alphas(trial));
    
  for csim = 1:nSim
    
    % Generate matrix based on spectrum
    fprintf("\n !!! Running simulation %d out of %d using alpha %d !!!\n", ...
      csim, nSim, alphas(trial));
    [Y, ~, Sigma] = synthetic_data_gen( feats, T, 1, alphas(trial) );
    
    % Test sapca edge
    if sapca_run == 1
        fprintf("\n\n -- Running SAPCA Edge test\n\n");

        fprintf("\n ** Running with rank seed of %d (out of: %d)\n", ...
          r_seed, feats);
        [~, fro_err_am(trial, csim, :), Uam, ~, ~, Yr, ...
          ... % for ticking
          t_sapca(csim, trial), ~, rmax, ...
          relerrs_am(trial,csim, :), rpcs_am(trial, csim, :)] ...
          = sapca_edge(Y, am_r_seed, am_tr_bounds, am_holdoff, am_blk, ...
          am_floor_mul, am_no_err, am_silent);

        % plot the fro error
        %plot(Tam, ErrFroAm, 'LineWidth', 2);
        am_yr = size(Yr, 2);
        %fprintf(" !! mse: %d (with seed %d)\n", T*mse(Yr, Y), rr_seed);
        fro_sapca_rseed = sum(sum((Y(:, 1:am_yr)-Yr).^2, 1))/am_yr;
        fprintf(" !! Final fro: %d (with seed %d)\n", ...
          fro_sapca_rseed, r_seed);

        errFroAMFinal(csim, trial) = fro_sapca_rseed;
        errMSEAMFinal(csim, trial) = am_yr*mse(Y(:, 1:am_yr), Yr);

        fprintf("\n\n -- Finished running SAPCA edge test\n");
    end

    % Test SPCA edge
    if spca_run == 1

      fprintf(" -- Running SPCA Edge test\n\n");

      lo_r = min(r_seed, size(Uam, 2));
      fprintf("\n ** Running low rank (r: %d)\n", lo_r);
      [~, fro_err_fm_lo(trial, csim, :), Ufm_lo, ~, ~, Yf_lo, ...
        ... % for ticking 
        t_spca_lo(csim, trial)] = ...
        spca_edge(Y, lo_r, fm_blk, fm_floor_mul, fm_no_err);

      lo_yr = size(Yf_lo, 2);
      fro_spca_lo_r = sum(sum((Y(:, 1:lo_yr)-Yf_lo).^2, 1))/lo_yr;
      fprintf(" !! Final fro: %d (low_r: %d)\n", fro_spca_lo_r, lo_r);

      errFroFMLoFinal(csim, trial) = fro_spca_lo_r;
      errMSEFMLoFinal(csim, trial) = lo_yr*mse(Y(:, 1:lo_yr), Yf_lo);

      % plot the low, all the time
      %plot(Tfm_lo, ErrFroFm_lo, '--');

      % need to run high rank to bound it
      if sapca_run == 1
        hi_r = rmax;
        fprintf("\n ** Running high rank (r: %d) to bound\n", hi_r);
        [~, fro_err_fm_hi(trial, csim, :), Ufm_hi, ~, ~, Yf_hi, ...
          ... % for ticking
          t_spca_hi(csim, trial)] = ...
          spca_edge(Y, hi_r, fm_blk, fm_floor_mul, fm_no_err);
        hi_yr = size(Yf_hi, 2);
        fro_spca_hi_r = sum(sum((Y(:, 1:hi_yr)-Yf_hi).^2, 1))/hi_yr;
        fprintf(" !! Final fro: %d (rank: %d)\n", fro_spca_hi_r, hi_r);

        %plot(Tfm_hi, ErrFroFm_hi);

        errFroFMHiFinal(csim, trial) = fro_spca_hi_r;
        errMSEFMHiFinal(csim, trial) = hi_yr*mse(Y(:, 1:hi_yr), Yf_hi);
      end

      fprintf("\n -- Finished running SPCA edge test\n");
    end

    % Test Mitliagkas Power Method
    if pm_run == 1
      pm_rank = min(pm_rank, lo_r);
      fprintf("\n -- Running PM with rank %d (out of: %d)\n", ...
        r_seed, feats);
      % ticking is embedded in the outputs
      [Tpm_lo, fro_err_pm_lo(trial, csim, :), Upms, YrPM, t_pm(csim, trial)] = ...
        mitliag_pm(Y, pm_rank, pm_blk, pm_floor_mul, pm_no_err);

      lo_yr = size(YrPM, 2);
      fro_pm_lo_r = sum(sum((Y(:, 1:lo_yr)-YrPM).^2, 1))/lo_yr;

      %plot(Tpm_lo, ErrFroPm_lo);
      % error calculation
      errFroPMFinal(csim, trial) = fro_pm_lo_r;
      errMSEPMFinal(csim, trial) = lo_yr*mse(Y(:, 1:lo_yr), YrPM);
      fprintf(" !! Final fro: %d (rank: %d)\n", fro_pm_lo_r, pm_rank);
      fprintf("\n -- Finished Running PM\n");
    end
    
    % Test FD
    if fd_run == 1
      fd_rank = min(fd_rank, lo_r);
      fprintf("\n -- Running FD with rank %d (out of %d)\n", ...
        r_seed, feats);
      % run fd
      [Ufd, fro_fd_lo(trial, csim, :), Tfd, Yr_fd, t_fd(csim, trial)] = ...
        fd(Y', fd_rank, fd_no_err);
      % transpose it
      Yr_fd = Yr_fd';
      Ufd = Ufd';
      
      lo_yr = size(Yr_fd, 2);
      fro_pm_lo_r = sum(sum((Y(:, 1:lo_yr)-Yr_fd).^2, 1))/lo_yr;
      
      errFroFDFinal(csim, trial) = fro_pm_lo_r;
      errMSEFDFinal(csim, trial) = lo_yr*mse(Y(:, 1:lo_yr), Yr_fd);
      fprintf(" !! Final fro: %d (rank: %d)\n", fro_pm_lo_r, fd_rank);
      fprintf("\n -- Finished Running FD\n");
    end
    
    % Test GROUSE
    if gr_run == 1
      gr_rank = min(gr_rank, lo_r);
      fprintf("\n -- Running GROUSE with rank %d (out of %d)\n", ...
        r_seed, feats);
      
      % run grouse
      [~, ~, U_gr, V_gr, t_gr(csim, trial)] = my_grouse(Y, gr_rank);
 
      % expand U_gr*V_gr' to get the Yr_gr
      Yr_gr = U_gr*V_gr';   
      gr_yr = size(Yr_gr, 2);
      % now get the fro error
      fro_gr_lo_r = sum(sum((Y(:, 1:gr_yr)-Yr_gr).^2, 1))/gr_yr;
      
      errFroGRFinal(csim, trial) = fro_gr_lo_r;
      errMSEGRFinal(csim, trial) = gr_yr*mse(Y(:, 1:gr_yr), Yr_fd);      
      fprintf(" !! Final fro: %d (rank %d)\n", fro_gr_lo_r, gr_rank);
      fprintf("\n -- Finished Running GROUSE\n");
    end


    % Test SPIRIT
    if sp_run == 1
      fprintf("\n -- Running SPIRIT test\n\n");
      
      % start the ticking
      sp_tic = tic;
      [W1,  k1, Proj1, recon1, ...
        relerrs_sp(trial, csim, :), rpcs_sp(trial, csim, :), ~, ...
        fro_err_sp(trial, csim, :)] = ...
        SPIRIT(Y', sp_lambda, sp_energy, sp_k0, sp_holdOffTime, ...
        sp_silent, sp_no_err);
      % finish the ticking
      t_spirit(csim, trial) = my_toc(sp_tic); 
      
      % error calculation
      YSpiritSubRecon = (W1*W1')*Y;
      fro_sp_rseed = sum(sum((Y-YSpiritSubRecon).^2, 1))/T;
      errFroSPFinal(csim, trial) = fro_sp_rseed;
      errMSESPFinal(csim, trial) = T*mse(Y, YSpiritSubRecon);
      fprintf(" !! Final fro: %d (rank: %d)\n", fro_sp_rseed, k1);
      fprintf("\n -- Finished Running SPIRIT\n");
    end

    %% Test the subspaces MSE

    if subspace_err_print == 1

      % Compute the offline Subspace of Y
      [Uoff, ~, ~] = svd(Y);
      % expand it, for comparison against the other subspaces
      pcaUU = Uoff*Uoff';
      % to have the absolute value
      Uoff_abs = abs(Uoff);
      % Subspace for SPIRIT
      if sp_run == 1
        %subTopRMSE = mse(pcaUU, W1*W1');
        subTopRMSE = mse(Uoff_abs(:, 1:k1), abs(W1));
        fprintf("\n ** SPIRIT Subspace (for r: %d) MSE: %d", ... 
          k1, subTopRMSE);
        subspaceTopRSPFinal(csim, trial) = subTopRMSE;
      end

      % Subspace for Power Method
      if pm_run == 1
        [~, r] = size(Upms);
        subTopRMSE = mse(Uoff_abs(:, 1:r), abs(Upms));
        fprintf("\n ** Mitliagkas Subspace (r: %d) MSE: %d", ...
          r, subTopRMSE);
        subspaceTopRPMFinal(csim, trial) = subTopRMSE;
      end
      
      % Subspace for FD
      if fd_run == 1
        [~, r] = size(Ufd);
        %subTopRMSE = mse(pcaUU, Ufd*Ufd');
        subTopRMSE = mse(Uoff_abs(:, 1:r), abs(Ufd));
        fprintf("\n ** FD Subspace (r: %d) MSE: %d", ...
          r, subTopRMSE);
        subspaceTopRFDFinal(csim, trial) = subTopRMSE;
      end
      
      % Subspace for GROUSE
      if gr_run == 1
        [~, r] = size(U_gr);
        subTopRMSE = mse(Uoff_abs(:, 1:r), abs(U_gr));
        fprintf("\n ** GROUSE Subspace (r: %d) MSE: %d", ...
          r, subTopRMSE);
        subspaceTopRGRFinal(csim, trial) = subTopRMSE;
      end

      % Subspace for SPCA Low
      if spca_run == 1
        % lower bound
        [~, r] = size(Ufm_lo);
        subTopRMSE = mse(Uoff_abs(:, 1:r), abs(Ufm_lo));
        %subTopRMSE = mse(pcaUU, Ufm_lo*Ufm_lo');
        fprintf("\n ** SPCA Subspace (lo_r: %d) MSE: %d", ...
          r, subTopRMSE);
        subspaceTopRFMLoFinal(csim, trial) = subTopRMSE;
        % higher bound
        [~, r] = size(Ufm_hi);
        subTopRMSE = mse(pcaUU, Ufm_hi*Ufm_hi');
        fprintf("\n ** SPCA Subspace UU' (hi_r: %d) MSE: %d", ...
          r, subTopRMSE);
        subspaceTopRFMHiFinal(csim, trial) = subTopRMSE;

      end

      % Subspace for SAPCA
      if sapca_run == 1
        [~, r] = size(Uam);
        subTopRMSE = mse(Uoff_abs(:, 1:r), abs(Uam));
        fprintf("\n ** SPCA Subspace (r: %d) MSE: %d", ...
          r, subTopRMSE);  
        subspaceTopRAMFinal(csim, trial) = subTopRMSE;
      end

    end
    % to end with a nice console offset
    fprintf("\n");
  
  end

end

% finally plot the aggregated errors

% print the Yr vs Y fro final error
if fro_print == 1
  fig = figure;
  hold on;
  plot(mean(errFroSPFinal, 1), '--');
  plot(mean(errFroAMFinal, 1), 'LineWidth', 2);
  plot(mean(errFroFMLoFinal, 1));
  plot(mean(errFroFMHiFinal, 1));
  plot(mean(errFroPMFinal, 1));
  plot(mean(errFroFDFinal, 1));
  hold off;
  if prod_print == 0
    stl = sprintf(['alphas vs fro errors ' , ...
      '(seed rank: %d, sims: %d, T: %d, feats: %d)'], ...
      r_seed, nSim, T, feats);
  else
    stl = sprintf('fro errors for varying alpha');
  end
  title(stl);
  xticks(1:alphas_len);
  xticklabels(alphas);
  ylabel('fro error');
  xlabel('alphas');
  legend('SP', 'SAPCA', 'SPCA_{lo}', 'SPCA_{hi}', ...
    'PM_{lo}', 'FD', 'Location', 'best');
  st = sprintf("fro_err_rseed_%d_sims_%d_T_%d", r_seed, nSim, T);
  print_fig(fig, st);
end

% print the Yr vs Y mse final error
if mse_print == 1
  fig = figure;
  hold on;
  plot(mean(errMSESPFinal, 1), '--');
  plot(mean(errMSEAMFinal, 1), 'LineWidth', 2);
  plot(mean(errMSEFMLoFinal, 1));
  plot(mean(errMSEFMHiFinal, 1));
  plot(mean(errMSEPMFinal, 1));
  plot(mean(errMSEFDFinal, 1));
  plot(mean(errMSEGRFinal, 1));
  hold off;
  stl = sprintf(['Yr vs real MSE ', ...
    '(seed rank: %d, sims: %d, T: %d, feats: %d)'], ...
    r_seed, nSim, T, feats);
  title(stl);
  xticks(1:alphas_len);
  xticklabels(alphas);
  ylabel('fro error');
  xlabel('alphas');
  legend('SP', 'SAPCA', 'SPCA_{lo}', 'SPCA_{hi}', ...
    'SPCA_{lo}', 'FD', 'GROUSE', 'Location', 'best');
  st = sprintf("mse_err_rseed_%d_sims_%d_T_%d", r_seed, nSim, T);
  print_fig(fig, st);
end

% print the subspace mse error
if subspace_err_print == 1
  fig = figure;
  hold on;
  plot(mean(subspaceTopRSPFinal, 1), '--');
  plot(mean(subspaceTopRAMFinal, 1), 'LineWidth', 2);
  plot(mean(subspaceTopRFMLoFinal, 1));
  plot(mean(subspaceTopRFMHiFinal, 1));
  plot(mean(subspaceTopRPMFinal, 1));
  plot(mean(subspaceTopRFDFinal, 1));
  plot(mean(subspaceTopRGRFinal, 1));
  hold off;
  if prod_print == 0
    stl = sprintf(['subspace vs real MSE ', ... 
      '(seed rank: %d, sims: %d, T: %d, feats: %d)'], ... 
      r_seed, nSim, T, feats);
  else
    stl = sprintf('Estimated vs Real Subspace MSE');
  end
  title(stl);
  xticks(1:alphas_len);
  xticklabels(alphas);
  ylabel('error (MSE)');
  xlabel('alphas');
  legend('SP', 'SAPCA', 'SPCA_{lo}', 'SPCA_{hi}', 'PM_{lo}', ...
    'FD', 'GROUSE', 'Location', 'SouthWest');
  % print to a file
  st = sprintf("subspace_err_feat_%d_T_%d_nsims_%d", ...
    feats, T, nSim);
  print_fig(fig, st);
end

% rel errors print
if rel_err_print == 1
  for i = 1:alphas_len
    fig = figure; 
    % pc evolution
    subplot(2, 1, 1);
      hold on;
      plot(mean(squeeze(rpcs_am(i, :, :)), 1), '--', 'LineWidth', 2); 
      plot(mean(squeeze(rpcs_sp(i, :, :)), 1), 'LineWidth', 2); 
      legend('AM', 'SP', 'Location', 'best'); 
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
      plot(mean(squeeze(relerrs_am(i, :, :)), 1), '--', 'LineWidth', 2);
      plot(mean(squeeze(relerrs_sp(i, :, :)), 1), 'LineWidth', 2);
      legend('SAPCA', 'SP', 'Location', 'best');
      stl = sprintf('Rel errors over T (%d) for alpha: %6.5f', T, alphas(i));
      title(stl);
      ylabel('rel. error');
      xlabel('time ticks');
      hold off;

      % print to a file
      st = sprintf("rel_errors_feat_%d_T_%d_alpha_%s_nsims_%d", ...
        feats, T, strrep(num2str(alphas(i)), ".", "_"), nSim);
    print_fig(fig, st);
  end
end

% rel errors print
if ext_fro_print == 1
  for i = 1:alphas_len
    fig = figure; 
    % pc evolution
    subplot(2, 1, 1);
      hold on;
      plot(mean(squeeze(rpcs_am(i, :, :)), 1), '--', 'LineWidth', 2); 
      plot(mean(squeeze(rpcs_sp(i, :, :)), 1), 'LineWidth', 2); 
      legend('SAPCA', 'SP', 'Location', 'best'); 
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
      plot(mean(squeeze(fro_err_am(i, :, :)), 1), '--', 'LineWidth', 2);
      plot(mean(squeeze(fro_err_sp(i, :, :)), 1), 'LineWidth', 2);
      legend('SAPCA', 'SP', 'Location', 'best');
      stl = sprintf('Fro errs over T (%d) for alpha: %6.5f', T, alphas(i));
      title(stl);
      ylabel('error (fro)');
      xlabel('time ticks');
      hold off;

      % print to a file
      st = sprintf("fro_ext_errors_feat_%d_T_%d_alpha_%s_nsims_%d", ...
        feats, T, strrep(num2str(alphas(i)), ".", "_"), nSim);
    print_fig(fig, st);
  end
end

% execution times print
if times_print == 1
  fig = figure;
  hold on;
  plot(mean(t_spirit, 1));
  plot(mean(t_sapca, 1));
  plot(mean(t_spca_lo, 1));
  plot(mean(t_spca_hi, 1));
  plot(mean(t_pm, 1));
  plot(mean(t_fd, 1));
  plot(mean(t_gr, 1));
  hold off;
  stl = sprintf("Mean Execution times (sims: %d, T: %d, feats: %d)", ...
    nSim, T, feats);
  title(stl);
  xticks(1:alphas_len);
  xticklabels(alphas);
  ylabel('execution time (s)');
  xlabel('alphas');
  % print to a file
  legend('SP', 'SAPCA', 'SPCA_{lo}', 'SPCA_{hi}', 'PM_{lo}', 'FD');
  st = sprintf("exec_times_feat_%d_T_%d_nsims_%d", ...
    feats, T, nSim);
  print_fig(fig, st);
end


