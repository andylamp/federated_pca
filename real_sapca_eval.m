function [err, legs] = real_sapca_eval(path, r, desc)
%REAL_SAPCA_EVAL function is responsible to run the comparison between 
% SAPCA, SPCA, GROUSE, SPIRIT, FD, and PM for a given real dataset pointed
% by `path` with starting rank `r`.
%
% Author: Andreas Grammenos (ag926@cl.cam.ac.uk)
%
% Last touched date: 18/07/2019
% 
% License: GPLv3
%

% scope in globals
global pflag
global fd_print

% initialisation
err = [];
legs = {};

% print iteration info
fprintf("\n\tTarget rank %d", r);
fprintf("\n\tPrint flag is: %d\n", pflag);

% load the dataset
Y = load(path)';

% perform alignment
block_pad = 50;    % round up to the nearest block that is multiple of this
[rows, cols] = size(Y);
bpad = mod(cols, block_pad);
Y = Y(:, 1: (cols-bpad));

% update the column number
cols = size(Y, 2);

% center the dataset by using Y = Y - ((Y*(vec_ones*vec_ones'))./cols)
vec_ones = ones(cols, 1);
Y = Y - ((Y*(vec_ones*vec_ones'))./cols);

% Parameters

% target rank (or seed for adaptive)
r_seed = r;
% enable analytical error calculation
no_err_flag = 0;

% SPIRIT parameters
sp_holdOffTime = 0;
sp_k0 = r_seed;
sp_lambda = .9;
sp_energy = [0.95, 0.98];
sp_silent = 1;
sp_no_err = no_err_flag;

% SAPCA parameters
am_r_seed = r_seed;
am_tr_bounds = [1, 10]; % percentages
am_holdoff = 0;
am_blk = 2*r_seed;
am_floor_mul = 2;
am_no_err = no_err_flag;
am_silent = 1;

% SPCA parameters
fm_rank = r_seed;
fm_blk = 2*r_seed;
fm_floor_mul = 2;
fm_no_err = no_err_flag;

% Frequent Directions Parameters
fd_rank = r_seed;
fd_err_flag = no_err_flag;

% Grouse Parameters
gr_rank = r_seed;
gr_no_err = no_err_flag;

% PM Parameters
pm_rank = r_seed;
pm_blk = rows;
pm_floor_mul = 2;
pm_no_err = no_err_flag;

% what to run
sapca_run = 1;
spca_run = 1;
pm_run = 1;
sp_run = 1;
fd_run = fd_print;
gr_run = 1;
grouse_run = 1;
fro_print = 1;
subspace_err_print = 1;

% production print, for shorter titles 
% used in the paper
prod_print = 1;

% Test SAPCA edge
if sapca_run == 1
    fprintf("\n\n -- Running SAPCA Edge test\n\n");

    fprintf("\n ** Running with rank seed of %d (out of: %d)\n", ...
      r_seed, rows);
    [Tam, ErrFroAm, Uam, ~, ~, Yr, ...
      ... % for ticking
      t_sapca, ~, rmax, ...
      relerrs_am, rpcs_am] ...
      = sapca_edge(Y, am_r_seed, am_tr_bounds, am_holdoff, am_blk, ...
      am_floor_mul, am_no_err, am_silent);

    am_yr = size(Yr, 2);
    fro_sapca_rseed = sum(sum((Y(:, 1:am_yr)-Yr).^2, 1))/am_yr;
    fprintf(" !! Final fro: %d (with seed %d)\n", ...
      fro_sapca_rseed, r_seed);

    fprintf("\n\n -- Finished running SAPCA edge test\n");
end

% Test SPCA edge
if spca_run == 1

  fprintf(" -- Running SPCA test\n\n");

  lo_r = min(fm_rank, size(Uam, 2));
  fprintf("\n ** Running low rank (r: %d)\n", lo_r);
  [Tfm_lo, ErrFroFm_lo, Ufm_lo, ~, ~, Yf_lo, ...
    ... % for ticking 
    t_spca_lo] = ...
    spca_edge(Y, lo_r, fm_blk, fm_floor_mul, fm_no_err);

  lo_yr = size(Yf_lo, 2);
  fro_spca_lo_r = sum(sum((Y(:, 1:lo_yr)-Yf_lo).^2, 1))/lo_yr;
  fprintf(" !! Final fro: %d (low_r: %d)\n", fro_spca_lo_r, lo_r);

  % need to run high rank to bound it
  if sapca_run == 1
    hi_r = rmax;
    fprintf("\n ** Running high rank (r: %d) to bound\n", hi_r);
    [Tfm_hi, ErrFroFm_hi, Ufm_hi, ~, ~, Yf_hi, ...
      ... % for ticking
      t_spca_hi] = ...
      spca_edge(Y, hi_r, fm_blk, fm_floor_mul, fm_no_err);
    hi_yr = size(Yf_hi, 2);
    fro_spca_hi_r = sum(sum((Y(:, 1:hi_yr)-Yf_hi).^2, 1))/hi_yr;
    fprintf(" !! Final fro: %d (rank: %d)\n", fro_spca_hi_r, hi_r);
  end

  fprintf("\n -- Finished running SPCA edge test\n");
end

% Test Mitliagkas Power Method
if pm_run == 1
  pm_rank = min(pm_rank, lo_r);
  fprintf("\n -- Running PM with rank %d (out of: %d)\n", ...
    pm_rank, rows);
  % ticking is embedded in the outputs
  [Tpm_lo, ErrFroPm_lo, Upms, YrPM, t_pm] = ...
    mitliag_pm(Y, pm_rank, pm_blk, pm_floor_mul, pm_no_err);

  lo_yr = size(YrPM, 2);
  fro_pm_lo_r = sum(sum((Y(:, 1:lo_yr)-YrPM).^2, 1))/lo_yr;

  fprintf(" !! Final fro: %d (rank: %d)\n", fro_pm_lo_r, pm_rank);
  fprintf("\n -- Finished Running PM\n");
end

% Test Frequent Directions
if fd_run == 1
  fprintf("\n -- Running FD with rank %d (out of: %d)\n", ...
    fd_rank, rows);
  % run the fd
  [Ufd, ErrFroFD, Tfd, Yr_fd, ~] = fd(Y', fd_rank, fd_err_flag);
  % since this is the transpose, revert it
  Yr_fd = Yr_fd';
  Ufd = Ufd';
  
  
  fd_yr = size(Yr_fd, 2);
  fro_fd_lo_r = sum(sum((Y(:, 1:fd_yr)-Yr_fd).^2, 1))/fd_yr;

  fprintf(" !! Final fro: %d (rank: %d)\n", fro_fd_lo_r, fd_rank);
  fprintf("\n -- Finished Running FD\n");
end

% Test Grouse
if grouse_run == 1
  fprintf("\n -- Running GROUSE with rank %d (out of: %d)\n", ...
    gr_rank, rows);
  
  [TGrouse, ErrFroGrouse, U_gr, V_gr, ~] = my_grouse(Y, gr_rank, gr_no_err);

  % expand U_gr*V_gr' to get the Yr_gr
  Yr_gr = U_gr*V_gr';
  
  grouse_yr = size(Yr_gr, 2);
  fro_grouse_lo_r = sum(sum((Y(:, 1:grouse_yr)-Yr_gr).^2, 1))/grouse_yr;
  
  fprintf(" !! Final fro: %d (rank: %d)\n", fro_grouse_lo_r, gr_rank);
  fprintf("\n -- Finished Running GROUSE\n");
  
end


% Test SPIRIT
if sp_run == 1
  fprintf("\n -- Running SPIRIT test\n\n");

  % start the ticking
  sp_tic = tic;
  [W1,  k1, Proj1, recon1, ...
    relerrs_sp, rpcs_sp, Tsp, ErrFroSP] = ...
    SPIRIT(Y', sp_lambda, sp_energy, sp_k0, sp_holdOffTime, sp_silent, ...
    sp_no_err);
  % finish the ticking
  t_spirit = my_toc(sp_tic); 

  % error calculation
  YSpiritSubRecon = (W1*W1')*Y;
  fro_sp_rseed = sum(sum((Y-YSpiritSubRecon).^2, 1))/cols;
  fprintf(" !! Final fro: %d (rank: %d)\n", fro_sp_rseed, k1);
  fprintf("\n -- Finished Running SPIRIT\n");
end

%% Test the subspaces MSE

if subspace_err_print == 1
  
  % Compute the offline PCA of Y
  [Upca, ~, ~] = svd(Y);
  % expand it, for comparison against the other subspaces
  %pcaUU = Upca*Upca';
  Uabs_off = abs(Upca);
  Usp = W1';
  idx = 1;
  % find the correct subspace rank to compare across all
  r = min([size(Upms, 2), k1, ...
    size(U_gr, 2), size(Ufm_lo, 2), size(Ufm_hi, 2)]);
  
  if fd_print == 1
    r = min([size(Ufd, 2), r]);
  end
  
  fprintf("\n Min r for all is %d\n", r); 
  
  
  % Subspace for SPIRIT
  
  if sp_run == 1
    %[~, r] = size(Usp);
    subspaceTopRSPFinal = mse(Uabs_off(:, 1:r), Usp);
    fprintf("\n ** SPIRIT Subspace (for r: %d) MSE: %d", ... 
      k1, subspaceTopRSPFinal);
    err(idx) = subspaceTopRSPFinal;
    legs{idx} = 'SP';
    idx = idx + 1;
  end 

  % Subspace for Power Method
  if pm_run == 1
    %[~, r] = size(Upms);
    subspaceTopRPMFinal = mse(Uabs_off(:, 1:r), Upms);
    fprintf("\n ** Mitliagkas Subspace (r: %d) MSE: %d", ...
      r, subspaceTopRPMFinal);
    err(idx) = subspaceTopRPMFinal;
    legs{idx} = 'PM';
    idx = idx + 1;
  end

  % Subspace for SPCA
  if spca_run == 1
    % lower bound
    %[~, r] = size(Ufm_lo);
    subspaceTopRFMLoFinal = mse(Uabs_off(:, 1:r), Ufm_lo);
    fprintf("\n ** SPCA Subspace (lo_r: %d) MSE: %d", ...
      r, subspaceTopRFMLoFinal);
    err(idx) = subspaceTopRFMLoFinal;
    legs{idx} = 'SPCA_{lo}';
    idx = idx + 1;
    % higher bound
    %[~, r] = size(Ufm_hi);
    subspaceTopRFMHiFinal = mse(Uabs_off(:, 1:r), Ufm_hi);
    fprintf("\n ** SPCA Subspace (hi_r: %d) MSE: %d", ...
      r, subspaceTopRFMHiFinal);
    err(idx) = subspaceTopRFMHiFinal;
    legs{idx} = 'SPCA_{hi}';
    idx = idx + 1;
  end

  % Subspace for FD
  if fd_run == 1
    %[~, r] = size(Ufd);
    subspaceTopRFDFinal = mse(Uabs_off(:, 1:r), Ufd);
    fprintf("\n ** FD Subspace (r: %d) MSE: %d", ...
      r, subspaceTopRFDFinal);
    err(idx) = subspaceTopRFDFinal;
    legs{idx} = 'FD';
    idx = idx + 1;
  end
  
  % Subspace for GROUSE
  if gr_run == 1
    %[~, r] = size(U_gr);
    subspaceTopRGRFinal = mse(Uabs_off(:, 1:r), U_gr);
    fprintf("\n ** GROUSE Subspace (r: %d) MSE: %d", ...
      r, subspaceTopRGRFinal);
    err(idx) = subspaceTopRGRFinal;
    legs{idx} = 'GROUSE';
    idx = idx + 1;
  end
  
  % Subspace for SPCA
  if sapca_run == 1
    %[~, r] = size(Uam);
    subspaceTopRAMFinal = mse(Uabs_off(:, 1:r), Uam);
    fprintf("\n ** SAPCA Subspace (r: %d) MSE: %d", ...
      r, subspaceTopRAMFinal);  
    err(idx) = subspaceTopRAMFinal;
    legs{idx} = 'SAPCA';
  end
  
end
% to end with a nice console offset
fprintf("\n");
  
% print the Yr vs Y fro final error
if fro_print == 1
  fig = figure;
  subplot(2, 1, 1);
  hold on;
  semilogy(Tsp, ErrFroSP, '--');
  semilogy(Tam, ErrFroAm, 'LineWidth', 2, ...
    'color', [0.9100    0.4100    0.1700]);
  semilogy(Tfm_lo, ErrFroFm_lo);
  semilogy(Tfm_hi, ErrFroFm_hi);
  semilogy(Tpm_lo, ErrFroPm_lo);
  semilogy(TGrouse, ErrFroGrouse);
  if fd_run == 1
    semilogy(Tfd, ErrFroFD);
  end
%     plot(ErrFroSP, Tsp, '--');
%     plot(ErrFroAm, Tam, 'LineWidth', 2);
%     plot(ErrFroFm_lo, Tfm_lo);
%     plot(ErrFroFm_hi, Tfm_hi);
%     plot(ErrFroPm_lo, Tpm_lo);
  hold off;
  if prod_print == 0
    stl = sprintf(['fro errors over time for %s Data' , ...
      '(seed rank: %d'], ...
      desc, r_seed);
  else
    stl = sprintf('fro errors for %s Data', desc);
  end
  title(stl);
  %xticks(1:alphas_len);
  %xticklabels(alphas);
  ylabel('fro error');
  xlabel('time ticks');
  legendCells ={'SP', 'SAPCA', 'SPCA_{lo}', 'SPCA_{hi}', 'PM_{lo}', ...
    'GROUSE', 'FD'};

  if fd_run == 0
    idc = ismember(legendCells, {'FD'});
    legendCells = legendCells(~idc);
  end
  legend(legendCells, 'Location', 'best');

  subplot(2, 1, 2);
  hold on;
  %semilogy(Tsp, ErrFroSP, '--');
  semilogy(Tam, ErrFroAm, 'LineWidth', 2, ...
    'color', [0.9100    0.4100    0.1700]);
  semilogy(Tfm_lo, ErrFroFm_lo);
  semilogy(Tfm_hi, ErrFroFm_hi);
  semilogy(Tpm_lo, ErrFroPm_lo);
%     plot(ErrFroSP, Tsp, '--');
%     plot(ErrFroAm, Tam, 'LineWidth', 2);
%     plot(ErrFroFm_lo, Tfm_lo);
%     plot(ErrFroFm_hi, Tfm_hi);
%     plot(ErrFroPm_lo, Tpm_lo);
  hold off;
  if prod_print == 0
    stl = sprintf(['fro errors over time for %s Data' , ...
      '(seed rank: %d'], ...
      desc, r_seed);
  else
    stl = sprintf('fro errors for %s Data without SP', desc);
  end
  title(stl);
  %xticks(1:alphas_len);
  %xticklabels(alphas);
  ylabel('fro error');
  xlabel('time ticks');
  legend('SAPCA', 'SPCA_{lo}', 'SPCA_{hi}', 'PM_{lo}', 'Location', 'best');
  %legend('AM', 'FM_{lo}', 'FM_{hi}', 'PM_{lo}', 'Location', 'best');
  st = sprintf("fro_err_rseed_real_%s", desc);
  % set the font size
  set(gca, 'FontSize', 12);
end
print_fig(fig, st);
end