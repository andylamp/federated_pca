function [err, legs] = eval_fpca_real(path, r, desc, params)
%EVAL_FPCA_REAL wrapper function to run all fpca real
%
% Based on work of Grammenos et al.: https://arxiv.org/abs/1907.08059
%
% Author: Andreas Grammenos (ag926@cl.cam.ac.uk)
%
% Last touched date: 13/06/2020
% 
% License: GPLv3
%

% check if we use block error instead (much faster)
if ~isfield(params, 'use_blk_err')
  params.use_blk_err = 0;
end

% check if we have printing enabled
if ~isfield(params, 'pflag')
  params.pflag = 0;
end

% check if we have a csv
if ~isfield(params, 'is_csv')
  params.is_csv = 0;
end

% check if we also compute frequent directions
if ~isfield(params, 'fd_run')
  params.fd_run = 0;
end

% check if we want to use a subplot or not
if ~isfield(params, 'use_subplot')
  params.use_subplot = 0;
end

% check if we have a production print mode or not
if ~isfield(params, 'prod_print')
  params.prod_print = 0;
end

% initialisation
err = [];
legs = {};

% print iteration info
fprintf("\n\tTarget rank %d", r);
fprintf("\n\tPrint flag is: %d\n", params.pflag);

% load the dataset, handle if normal or csv
if params.is_csv == 1
  Y = csvread(path, 1)';
else
  Y = load(path)';
end

% perform alignment
block_pad = 50;    % round up to the nearest block that is multiple of this
[rows, cols] = size(Y);
bpad = mod(cols, block_pad);
Y = Y(:, 1: (cols-bpad));

% update the column number
cols = size(Y, 2);

% center & normalise by using Y = Y - ((Y*(vec_ones*vec_ones'))./cols)
vec_ones = ones(cols, 1);
Y = Y - ((Y*(vec_ones*vec_ones'))./cols);

% Parameters

% target rank (or seed for adaptive)
r_seed = r;
% enable analytical error calculation
no_err_flag = 0;


% SPIRIT parameters
sp_lambda = .9;
sp_energy = [.95, .98];

sp_params.verbose = 0;
sp_params.no_err = no_err_flag;
sp_params.use_blk_err = params.use_blk_err;
sp_params.k0 = r_seed;
sp_params.holdoff_time = 0;

% parameters for fpca

% adaptive
fpca_params.adaptive = 1;
fpca_params.blk_size = 2*r_seed;
fpca_params.no_err = no_err_flag;
fpca_params.use_blk_err = params.use_blk_err;

% fixed (low rank)
fpca_low_params.adaptive = 0;
fpca_low_params.blk_size = 2*r_seed;
fpca_low_params.no_err = no_err_flag;
fpca_low_params.use_blk_err = params.use_blk_err;

% fixed (high rank)
fpca_high_params.adaptive = 0;
fpca_high_params.blk_size = 2*r_seed;
fpca_high_params.no_err = no_err_flag;
fpca_high_params.use_blk_err = params.use_blk_err;

% Frequent Directions Parameters
fd_rank = r_seed;
fd_params.use_blk_err = params.use_blk_err;
fd_params.no_err = no_err_flag;

% Grouse Parameters
gr_rank = r_seed;
gr_params.no_err = no_err_flag;
gr_params.use_blk_err = params.use_blk_err;

% PM Parameters
pm_rank = r_seed;
pm_params.no_err = no_err_flag;
pm_params.blk_size = rows;

% Test F-PCA (adaptive)
fprintf("\n\n -- Running F-PCA (adaptive) test\n\n");

% compute the adaptive fpca
[Ufpca, ~, fpca_opt_out] = fpca_edge(Y, r_seed, fpca_params);

yr_sz = size(fpca_opt_out.Yr, 2);
% fprintf(" !! mse: %d (with seed %d)\n", T*mse(Yr, Y), rr_seed);
fpca_fro = sum(sum((Y(:, 1:yr_sz)-fpca_opt_out.Yr).^2, 1))/ yr_sz;
fprintf(" !! Final fro: %d (with seed %d)\n", fpca_fro, r_seed);

fprintf("\n\n -- Finished running F-PCA (adaptive) test\n");

% Test F-PCA with adaptive behaviour disabled
if params.fpca_fixed_run == 1
  fprintf("\n\n -- Running F-PCA (fixed) test\n\n");
  
  % get the low rank, which is the min of the initial seed and the
  % recovered rank.
  lo_r = min(r_seed, size(Ufpca, 2));
  
  fprintf("\n ** Running low rank (r: %d)\n", lo_r);
  
  % compute the fpca for the min observed rank
  [Ufpca_lo, ~, fpca_low_opt_out] = fpca_edge(Y, lo_r, fpca_low_params);

  % compute the frobenius error
  fpca_lo_yr = size(fpca_low_opt_out.Yr, 2);
  yr_low_csum = sum(sum((Y(:, 1:fpca_lo_yr) - fpca_low_opt_out.Yr).^2, 1));
  fro_spca_lo_r =  yr_low_csum / fpca_lo_yr;
  
  fprintf(" !! Final fro: %d (low_r: %d)\n", fro_spca_lo_r, lo_r);

  % grab the highest rank observed
  hi_r = fpca_opt_out.rmax;

  fprintf("\n ** Running high rank (r: %d) to bound\n", hi_r);

  % compute the fpca for the max observed rank
  [Ufpca_hi, ~, fpca_high_opt_out] = fpca_edge(Y, hi_r, fpca_high_params);

  % compute the frobenious error
  hi_yr = size(fpca_high_opt_out.Yr, 2);
  yr_high_csum = sum(sum((Y(:, 1:hi_yr)-fpca_high_opt_out.Yr).^2, 1));
  fro_spca_hi_r =  yr_high_csum / hi_yr;
  fprintf(" !! Final fro: %d (rank: %d)\n", fro_spca_hi_r, hi_r);

  fprintf("\n -- Finished running F-PCA Fixed edge test\n");
end

% Test Mitliagkas Power Method
if params.pm_run == 1
  % pm_rank = min(pm_rank, r_seed);
  
  fprintf("\n -- Running PM with rank %d (out of: %d)\n", ...
    pm_rank, rows);
  
  % compute the power method
  [Upms, pm_opt_out] = mitliag_pm(Y, r_seed, pm_params);
  
  % compute the frobenious error for power method
  pm_yr_sz = size(pm_opt_out.Yr, 2);
  fro_pm_lo_r = sum(sum((Y(:, 1:pm_yr_sz)-pm_opt_out.Yr).^2, 1)) / pm_yr_sz;
  fprintf(" !! Final fro: %d (rank: %d)\n", fro_pm_lo_r, pm_rank);
  
  fprintf("\n -- Finished Running PM\n");
end

% Test Frequent Directions (only if enabled as it skews plots)
if params.fd_run == 1
  fprintf("\n -- Running FD with rank %d (out of: %d)\n", ...
    fd_rank, rows);
  
  % run frequent directions  
  [Ufd, fd_opt_out] = fd(Y', fd_rank, fd_params);
  
  % since this is the transpose, revert it
  Yr_fd = fd_opt_out.Yr';
  Ufd = Ufd';
  
  % compute the frobenious error
  fd_yr = size(Yr_fd, 2);
  fro_fd_lo_r = sum(sum((Y(:, 1:fd_yr)-Yr_fd).^2, 1))/fd_yr;

  fprintf(" !! Final fro: %d (rank: %d)\n", fro_fd_lo_r, fd_rank);
  fprintf("\n -- Finished Running FD\n");
end

% Test Grouse
if params.gr_run == 1
  fprintf("\n -- Running GROUSE with rank %d (out of: %d)\n", ...
    gr_rank, rows);
  
  % compute grouse
  [U_gr, V_gr, gr_opt_out] = my_grouse(Y, gr_rank, gr_params);

  % expand U_gr*V_gr' to get the Yr_gr
  Yr_gr = U_gr*V_gr';
  
  % compute the frobenious error
  grouse_yr = size(Yr_gr, 2);
  fro_grouse = sum(sum((Y(:, 1:grouse_yr)-Yr_gr).^2, 1)) / grouse_yr;

  fprintf(" !! Final fro: %d (rank: %d)\n", fro_grouse, gr_rank);
  fprintf("\n -- Finished Running GROUSE\n");
  
end

% Test SPIRIT
if params.sp_run == 1
  fprintf("\n -- Running SPIRIT test\n\n"); 

  % compute spirit
  [Usp, sp_opt_out] = SPIRIT(Y', sp_lambda, sp_energy, sp_params);

  % frobenious error calculation
  YSpiritSubRecon = (Usp*Usp')*Y;
  fro_sp_rseed = sum(sum((Y-YSpiritSubRecon).^2, 1))/cols;
  
  fprintf(" !! Final fro: %d (rank: %d)\n", ...
    fro_sp_rseed, sp_opt_out.final_rank);
  fprintf("\n -- Finished Running SPIRIT\n");
end

%% Test the subspaces MSE

if params.subspace_err_print == 1
  
  % Compute the offline PCA of Y
  [Upca, ~, ~] = svd(Y);
  % expand it, for comparison against the other subspaces
  %pcaUU = Upca*Upca';
  Uabs_off = abs(Upca);
  % Usp = W1';
  idx = 1;
  
  % get the subspace rank
  s_r = r;
  
  
  % find the correct subspace rank to compare based on execution parameters
  if params.sp_run == 1
    s_r = min([sp_opt_out.final_rank, s_r]);
  end
  
  % check the min for fpca
  if params.fpca_fixed_run == 1
    s_r = min([size(Ufpca_lo, 2), s_r]);
  else
    s_r = min([size(Ufpca, 2), s_r]);
  end
  
  
  fprintf("\n == Min rank (r) for all methods is %d\n", r); 
  
  % Subspace for SPIRIT
  if params.sp_run == 1
    subspaceTopRSPFinal = immse(Uabs_off(:, 1:s_r), Usp(:, 1:s_r));
    fprintf("\n ** SPIRIT Subspace (for r: %d) MSE: %d", ... 
      sp_opt_out.final_rank, subspaceTopRSPFinal);
    err(idx) = subspaceTopRSPFinal;
    legs{idx} = 'SP';
    idx = idx + 1;
  end 

  % Subspace for Power Method
  if params.pm_run == 1
    subspaceTopRPMFinal = immse(Uabs_off(:, 1:s_r), abs(Upms(:, 1:s_r)));
    fprintf("\n ** Mitliagkas Subspace (r: %d) MSE: %d", ...
      s_r, subspaceTopRPMFinal);
    err(idx) = subspaceTopRPMFinal;
    legs{idx} = 'PM';
    idx = idx + 1;
  end
  
  % Subspace for FD
  if params.fd_run == 1
    subspaceTopRFDFinal = immse(Uabs_off(:, 1:s_r), abs(Ufd(:, 1:s_r)));
    fprintf("\n ** FD Subspace (r: %d) MSE: %d", ...
      s_r, subspaceTopRFDFinal);
    err(idx) = subspaceTopRFDFinal;
    legs{idx} = 'FD';
    idx = idx + 1;
  end
  
  % Subspace for GROUSE
  if params.gr_run == 1
    subspaceTopRGRFinal = immse(Uabs_off(:, 1:s_r), abs(U_gr(:, 1:s_r)));
    fprintf("\n ** GROUSE Subspace (r: %d) MSE: %d", ...
      s_r, subspaceTopRGRFinal);
    err(idx) = subspaceTopRGRFinal;
    legs{idx} = 'GROUSE';
    idx = idx + 1;
  end
  
  % Subspace for F-PCA (adaptive)
  subspaceTopRAMFinal = immse(Uabs_off(:, 1:s_r), abs(Ufpca(:, 1:s_r)));
  fprintf("\n ** F-PCA Subspace (r: %d) MSE: %d", ...
    s_r, subspaceTopRAMFinal);  
  err(idx) = subspaceTopRAMFinal;
  legs{idx} = 'F-PCA';

  idx = idx + 1;
  
  % Subspace for F-PCA (fixed)
  if params.fpca_fixed_run == 1
    % lower bound
    subspaceTopRFMLoFinal = immse(Uabs_off(:, 1:s_r), ...
      abs(Ufpca_lo(:, 1:s_r)));
    fprintf("\n ** F-PCA (fixed) Subspace (lo_r: %d) MSE: %d", ...
      s_r, subspaceTopRFMLoFinal);
    err(idx) = subspaceTopRFMLoFinal;
    legs{idx} = 'F-PCA_{lo}';
    idx = idx + 1;
    % higher bound
    subspaceTopRFMHiFinal = immse(Uabs_off(:, 1:s_r), ...
      abs(Ufpca_hi(:, 1:s_r)));
    fprintf("\n ** F-PCA (fixed) Subspace (hi_r: %d) MSE: %d", ...
      s_r, subspaceTopRFMHiFinal);
    err(idx) = subspaceTopRFMHiFinal;
    legs{idx} = 'F-PCA_{hi}';
  end
  
end
% to end with a nice console offset
fprintf("\n");

% check if we indeed print
if params.fro_print == 1
  fig = figure;

  hold on;
  if params.sp_run == 1
    semilogy(sp_opt_out.T, log(sp_opt_out.ErrFro),'--');
  end

  if params.pm_run == 1
    semilogy(pm_opt_out.T, log(pm_opt_out.ErrFro));
  end

  if params.gr_run == 1
    semilogy(gr_opt_out.T, log(gr_opt_out.ErrFro));
  end

  if params.fd_run == 1
    semilogy(fd_opt_out.T, log(fd_opt_out.ErrFro));
  end

  semilogy(fpca_opt_out.T, log(fpca_opt_out.ErrFro), 'LineWidth', 2, ...
    'color', [0.4660    0.6740    0.1880]);

  if params.fpca_fixed_run == 1
    semilogy(fpca_low_opt_out.T, log(fpca_low_opt_out.ErrFro));
    semilogy(fpca_high_opt_out.T, log(fpca_high_opt_out.ErrFro));
  end
  hold off;

  % check if we have production print
  if params.prod_print == 0
    stl = sprintf(['fro errors over time for %s Data' , ...
      '(seed rank: %d'], ...
      desc, r_seed);
  else
    stl = sprintf('fro errors for %s Data', desc);
  end
  title(stl);

  ylabel('error (log(fro))');
  xlabel('time ticks');

  legendCells ={'SP', 'PM', 'GROUSE', 'FD', ...
    'F-PCA', 'F-PCA_{lo}', 'F-PCA_{hi}', };

  if params.sp_run == 0
    idc = ismember(legendCells, {'SP'});
    legendCells = legendCells(~idc);
  end

  if params.fd_run == 0
    idc = ismember(legendCells, {'FD'});
    legendCells = legendCells(~idc);
  end

  if params.pm_run == 0
    idc = ismember(legendCells, {'PM'});
    legendCells = legendCells(~idc);
  end

  if params.gr_run == 0
    idc = ismember(legendCells, {'GROUSE'});
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

  legend(legendCells, 'Location', 'best');

  % set the font size
  set(gca, 'FontSize', 12);

  % the title
  st = sprintf("fro_err_rseed_real_%s", desc);   

  % print the figure, if needed
  print_fig(fig, st, params);

end