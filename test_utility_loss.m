% The script is designed to evaluate the utility loss between F-PCA with
% perturbation masks against MOD-SuLQ (symmetric & non-symmetric).
%
% In all instances adaptive rank estimation is disabled for F-PCA.
%
% Based on work of Grammenos et al.: https://arxiv.org/abs/1907.08059
%
% Author: Andreas Grammenos (ag926@cl.cam.ac.uk)
%
% Last touched date: 31/05/2020
% 
% License: GPLv3
%

%% Initialisation

clc; clear; close all;

% for reproducibility
rng(300, 'twister');

% setup the variables
ul_params.type = "utility-loss";
% enable printing
ul_params.pflag = 1;
% print pdfs
ul_params.pdf_print = 1;
% configure the parameters
ul_params = setup_vars(ul_params);

% number of simulations to run
sims = 100;

% target rank (or seed for adaptive)
r = 5;
%rr = 2;
rr = 1;
% number of features in each vector
feats = 20; % 1k, 200 - for paper we used 20
% number of feature vectors to process
T = 5000; % 10k, 2k - for paper we used 5k
% singular value scaler
sv_scaler = 1;
% the alpha value range
alphas = [0.01, 0.1, .5, 1];

% show spectrums
show_spectra = 1;

% check if we use normal distribution
use_normal = 0;

% the privacy deltas
deltas = [.05]; % .01 
% the epsilon values
epsilons = 0.1:0.1:4; % 3

synth_params.spectrum_type = "pl";
%synth_params.spectrum_type = "singular";
%synth_params.spectrum_type = "mod-sulq";

% F-PCA parameters 
adaptive = 0;
private = 1;

% loop over the alphas

a_len = size(alphas, 2);
e_len = size(epsilons, 2);
d_len = size(deltas, 2);

% preallocate trade-off arrays
u_tradeoff = zeros(a_len, d_len, e_len);
ms_u_tradeoff = zeros(a_len, d_len, e_len);
ms_nonsym_u_tradeoff = zeros(a_len, d_len, e_len);

% pre-allocate the spectrum arrays only if we plot them
if show_spectra == 1
  sv_spectra = zeros(a_len, feats);
end

% preallocate cells for legends
ll = cell(a_len, 1); 

%% Run the utility loss experiments

for j = 1:size(alphas, 2)
  % configure power law, if it is the elected distribution type.
  if synth_params.spectrum_type == "pl"
    synth_params.alpha = alphas(j);
    synth_params.lambda = sv_scaler;
    synth_params.rand_type = "normal";
  end
  
  % generate the synthetic data
  if show_spectra == 1
    [Y, sv_spectra(j, :)] = synthetic_data_gen(feats, T, synth_params);
  else
    [Y, ~] = synthetic_data_gen(feats, T, synth_params);
  end
  [UamY, ~, ~] = svds(Y, rr);
  
  % get the correct amount of PC's to compare against (normally the first)
  Uam = UamY(:, 1:rr);
  
  % compute the Vperp
  Vperp = null(Uam');     
  
  % run the for number of simulations
  for kk = 1:sims
    % run for the delta range
    for d = 1:size(deltas, 2)
      delta = deltas(d);
      % run for the epsilon range
      for i = 1:size(epsilons, 2)
        e_p = epsilons(i);

        % normal mod-sulq
        [Ums, ~, ~] = svds(mod_sulq(Y, e_p, delta), rr);
        Ums = Ums(:, 1:rr);
        
        % non-symmetric mod-sulq
        [Ums_nonsym, ~, ~] = svds(mod_sulq(Y, e_p, delta, 0), rr);
        Ums_nonsym = Ums_nonsym(:, 1:rr);
        
         
        % configure f-pca
        params.adaptive = adaptive;
        params.private = private;
        params.e_p = e_p; 
        params.delta = delta;
        params.blk_size = 50;
        
        % run f-pca
        [Uam_p, ~, ~] = fpca_edge(Y, r, params);         
        Uam_p = Uam_p(:, 1:rr);        
        
        % compute the trade off for f-pca
        u_tradeoff(j, d, i) = u_tradeoff(j, d, i) + norm(Uam_p'*UamY);

        % non symmetric mod-sulq
        ms_nonsym_u_tradeoff(j, d, i) = ms_nonsym_u_tradeoff(j, d, i) + norm(Ums_nonsym'*UamY);
        
        % compute the trade-off for mod-sulq
        ms_u_tradeoff(j, d, i) = ms_u_tradeoff(j, d, i) + norm(Ums'*UamY);
      end
    end
  end
end


%%  Finally, plot the results
  
fig = figure;

% adjust the number of figures based if we show spectra or not
if show_spectra == 1
  plot_num = 4;
else
  plot_num = 3;
end

subplot(plot_num, 1, 1)
hold on
% run through the alphas and plot the results
for i = 1:a_len
  ff = u_tradeoff(i, :, :);
  plot(squeeze(ff)/sims, 'LineWidth', 2);
  s = sprintf("d at a %.2f", alphas(i));
  ll{i} = sprintf('a=%.2f', alphas(i));
end
hold off
% this is to output the "\varepsilon" equivalent
xlabel(char(949));
ylabel('q_{A}');

% configure the xticks
xlabels = xticks;
xlabels_sz = size(xlabels, 2);
xtick_labels = cell(1, xlabels_sz);
xtick_labels{1} = epsilons(1);
% assign the rest of the tick labels
for i = 2:xlabels_sz
  xtick_labels{i} = epsilons(xlabels(i));
end

% assign the ticks
xticklabels(xtick_labels);

title('fpca');
legend(ll);


subplot(plot_num, 1, 2)
hold on
% run through the alphas and plot the results
for i = 1:a_len
  ff = ms_nonsym_u_tradeoff(i, :, :);
  plot(squeeze(ff)/sims, 'LineWidth', 2);
  s = sprintf("d at a %.2f", alphas(i));
  ll{i} = sprintf('a=%.2f', alphas(i));
end
hold off
% this is to output the "\varepsilon" equivalent
xlabel(char(949));
ylabel('q_{A}');

% assign the ticks
xticklabels(xtick_labels);

title('mod-sulq (non-symmetric)');  
legend(ll);

subplot(plot_num, 1, 3)
hold on
% run through the alphas and plot the results
for i = 1:a_len
  ff = ms_u_tradeoff(i, :, :);
  plot(squeeze(ff)/sims, 'LineWidth', 2);
  s = sprintf("d at a %.2f", alphas(i));
  ll{i} = sprintf('a=%.2f', alphas(i));
end
hold off
% this is to output the "\varepsilon" equivalent
xlabel(char(949));
ylabel('q_{A}');

% assign the ticks
xticklabels(xtick_labels);

title('mod-sulq (symmetric)');  
legend(ll);


% check if we plot spectra and plot if we enabled.
if show_spectra == 1
  subplot(plot_num, 1, 4)
  hold on
  for i = 1:a_len
    ff = squeeze(sv_spectra(i, :));
    plot(ff, 'LineWidth', 2);
  end
  hold off
  xlabel('singular values')
  ylabel('magnitude');
  title('Singular Value True Spectrum')
  legend(ll);
end

% configure filename for the figure
st = sprintf("utility_loss_spectrum_type_%s_d_%d_r_%d_T_%d_sims_%d", ...
  synth_params.spectrum_type, feats, r, T, sims);
% print the figure for the utility loss
print_fig(fig, st, ul_params);


%% Plot them individually

fig = figure;
hold on
% run through the alphas and plot the results
for i = 1:a_len
  ff = u_tradeoff(i, :, :);
  plot(squeeze(ff)/sims, 'LineWidth', 2);
  s = sprintf("d at a %.2f", alphas(i));
  ll{i} = sprintf('a=%.2f', alphas(i));
end
hold off
% this is to output the "\varepsilon" equivalent
xlabel(char(949));
ylabel('q_{A}');

% configure the xticks
xlabels = xticks;
xlabels_sz = size(xlabels, 2);
xtick_labels = cell(1, xlabels_sz);
xtick_labels{1} = epsilons(1);
% assign the rest of the tick labels
for i = 2:xlabels_sz
  xtick_labels{i} = epsilons(xlabels(i));
end

% assign the ticks
xticklabels(xtick_labels);

title('fpca');
legend(ll, 'location', 'best');

% configure filename for the figure
st = sprintf("fpca_utility_loss_paper_spectrum_type_%s_d_%d_r_%d_T_%d_sims_%d", ...
  synth_params.spectrum_type, feats, r, T, sims);
% print the figure for the utility loss
print_fig(fig, st, ul_params);


fig = figure;
hold on
% run through the alphas and plot the results
for i = 1:a_len
  ff = ms_nonsym_u_tradeoff(i, :, :);
  plot(squeeze(ff)/sims, 'LineWidth', 2);
  s = sprintf("d at a %.2f", alphas(i));
  ll{i} = sprintf('a=%.2f', alphas(i));
end
hold off
% this is to output the "\varepsilon" equivalent
xlabel(char(949));
ylabel('q_{A}');

% assign the ticks
xticklabels(xtick_labels);

title('mod-sulq (non-symmetric)');  
legend(ll, 'location', 'best');

% configure filename for the figure
st = sprintf("mod_sulq_sym_utility_loss_paper_spectrum_type_%s_d_%d_r_%d_T_%d_sims_%d", ...
  synth_params.spectrum_type, feats, r, T, sims);
% print the figure for the utility loss
print_fig(fig, st, ul_params);

fig = figure;
hold on
% run through the alphas and plot the results
for i = 1:a_len
  ff = ms_u_tradeoff(i, :, :);
  plot(squeeze(ff)/sims, 'LineWidth', 2);
  s = sprintf("d at a %.2f", alphas(i));
  ll{i} = sprintf('a=%.2f', alphas(i));
end
hold off
% this is to output the "\varepsilon" equivalent
xlabel(char(949));
ylabel('q_{A}');

% assign the ticks
xticklabels(xtick_labels);

title('mod-sulq (symmetric)');  
legend(ll, 'location', 'best');

% configure filename for the figure
st = sprintf("mod_sulq_utility_loss_paper_spectrum_type_%s_d_%d_r_%d_T_%d_sims_%d", ...
  synth_params.spectrum_type, feats, r, T, sims);
% print the figure for the utility loss
print_fig(fig, st, ul_params);

