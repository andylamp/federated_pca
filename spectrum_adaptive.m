% Plot the estimated singular vs the ground truth ones. 
% Values were taken at dry runs and paste in this file for easier plotting.

% initialise
clc; clear; close all;

% for alpha 0.5
fpca_sv_05 =	[1.0069 0.7078 0.5717 0.5066 0.4451 ...
               0.4044 0.3801 0.3526 0.3276 0.3118];
	
real_sv_05 = 	[1.0000 0.7071 0.5774 0.5000 0.4472 ...
               0.4082 0.3780 0.3536 0.3333 0.3162];

% for alpha 1
fpca_sv_1 = [1.0067 0.5004 0.3300 0.2535 0.1992 ...
              0.1653 0.1436 0.1258 0.1105 0.0964];

real_sv_1 = [1.0000 0.5000 0.3333 0.2500 0.2000 ...
             0.1667 0.1429 0.1250 0.1111 0.1000];

% for alpha 2
fpca_sv_2 = [1.0067 0.2502 0.1100 0.0634 0.0398 ...
              0.0275 0.0205 0.0157 0.0014 0.0010];
            
real_sv_2 = [1.0000 0.2500 0.1111 0.0625 0.0400 ...
             0.0278 0.0204 0.0156 0.0123 0.0100];
           
% plot it nicely
r = size(fpca_sv_05, 2);
figure;
hold on;
plot(fpca_sv_05, '-*', 'LineWidth', 2);
plot(real_sv_05, '-+');
plot(fpca_sv_1, '-*', 'LineWidth', 2);
plot(real_sv_1, '-+');
plot(fpca_sv_2, '-*', 'LineWidth', 2);
plot(real_sv_2, '-+');
hold off;
title('Adaptive Singular value significance');
legend('F-PCA_{a=.5}', 'GT_{a=.5}', 'F-PCA_{a=1}', ...
  'GT_{a=1}', 'F-PCA_{a=2}', 'GT_{a=2}');
xlabel('Singular values');
xticklabels(1:r);
ylabel('Values');
