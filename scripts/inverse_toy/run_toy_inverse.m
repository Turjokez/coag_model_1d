% run_toy_inverse.m
%
% Toy inverse problem: fit 2 parameters to noisy depth-profile data.
%
% Model: P(z) = P0 * exp(-b * z)
%   P0 = surface concentration
%   b  = attenuation coefficient (1/m)
%
% Method:
%   1. Generate synthetic "observations" from true parameters + noise
%   2. Set prior values (different from true, to show recovery)
%   3. Minimize cost: J = data misfit + prior penalty
%   4. Estimate uncertainty from approximate Hessian
%
% Steps:
%   (a) cost function on a grid (see where the minimum is)
%   (b) fminsearch to find best fit
%   (c) finite-difference Hessian for posterior uncertainty
%   (d) two figures: cost surface, and data vs fit

set(0,'DefaultAxesFontName','Arial');
set(0,'DefaultTextFontName','Arial');

% ---------------------------------------------------------------
% 1. True parameters and synthetic observations
% ---------------------------------------------------------------

P0_true = 5.0;    % true surface value  (arbitrary units)
b_true  = 0.008;  % true attenuation    (m^-1)

z_obs   = [25, 75, 150, 250, 350, 500]';   % observation depths (m)
noise   = 0.15;                             % relative noise level

rng(42);
P_true = P0_true * exp(-b_true * z_obs);
P_obs  = P_true .* (1 + noise * randn(size(z_obs)));
P_obs  = max(P_obs, 0);    % keep positive

% observation error (std), same as noise level * true value
sigma_d = noise * P_true;

% ---------------------------------------------------------------
% 2. Prior values and uncertainties
% ---------------------------------------------------------------

% prior: deliberately offset from truth
P0_prior = 3.0;
b_prior  = 0.005;

% prior uncertainty (1 sigma)
sigma_P0 = 2.0;    % large -> prior is weak
sigma_b  = 0.004;

% ---------------------------------------------------------------
% 3. Cost function
% ---------------------------------------------------------------
% J = sum( (P_obs - P_mod)^2 / sigma_d^2 )
%   + (P0 - P0_prior)^2 / sigma_P0^2
%   + (b  - b_prior)^2  / sigma_b^2

cost = @(p) cost_fn(p, z_obs, P_obs, sigma_d, ...
                    P0_prior, b_prior, sigma_P0, sigma_b);

% ---------------------------------------------------------------
% 4. Grid search: see the cost surface
% ---------------------------------------------------------------

P0_vec = linspace(1, 10, 60);
b_vec  = linspace(0.002, 0.016, 60);
[P0g, bg] = meshgrid(P0_vec, b_vec);
Jg = zeros(size(P0g));

for i = 1:numel(P0g)
    Jg(i) = cost([P0g(i), bg(i)]);
end

% ---------------------------------------------------------------
% 5. fminsearch from prior as starting point
% ---------------------------------------------------------------

p0   = [P0_prior, b_prior];
opts = optimset('TolX', 1e-8, 'TolFun', 1e-10, 'MaxFunEvals', 2000);
[p_fit, J_fit] = fminsearch(cost, p0, opts);

P0_fit = p_fit(1);
b_fit  = p_fit(2);

fprintf('True:   P0 = %.3f,  b = %.5f\n', P0_true, b_true);
fprintf('Prior:  P0 = %.3f,  b = %.5f\n', P0_prior, b_prior);
fprintf('Fit:    P0 = %.3f,  b = %.5f   (J = %.4f)\n', P0_fit, b_fit, J_fit);

% ---------------------------------------------------------------
% 6. Approximate uncertainty from finite-difference Hessian
%    H ~ Jacobian^T Cd^-1 Jacobian + Cp^-1
%    C_post ~ inv(H)
%    Only data term Jacobian here (simpler).
% ---------------------------------------------------------------

dp    = [1e-4 * P0_fit, 1e-6 * b_fit];    % small steps
J_mat = zeros(numel(z_obs), 2);            % Jacobian of model at p_fit

for j = 1:2
    pp = p_fit; pp(j) = pp(j) + dp(j);
    pm = p_fit; pm(j) = pm(j) - dp(j);
    J_mat(:,j) = (fwd(pp, z_obs) - fwd(pm, z_obs)) / (2*dp(j));
end

Cd_inv = diag(1 ./ sigma_d.^2);
Cp_inv = diag([1/sigma_P0^2, 1/sigma_b^2]);

H      = J_mat' * Cd_inv * J_mat + Cp_inv;
C_post = inv(H);
sigma_post = sqrt(diag(C_post));

fprintf('Post sigma:  P0 = %.3f,  b = %.5f\n', sigma_post(1), sigma_post(2));

% ---------------------------------------------------------------
% 7. Figure 1: cost surface
% ---------------------------------------------------------------

fs = 7;
figure('Units','centimeters','Position',[2 2 10 8],'Color','white');
contour(P0g, bg*1000, log10(Jg), 20, 'LineWidth', 0.8, 'HandleVisibility', 'off');
hold on;
plot(P0_prior,  b_prior*1000,  'bs', 'MarkerSize', 5, 'DisplayName', 'prior');
plot(P0_fit,    b_fit*1000,    'r^', 'MarkerSize', 5, 'DisplayName', 'fit');
plot(P0_true,   b_true*1000,   'ko', 'MarkerSize', 5, 'DisplayName', 'true');
xlabel('P_0', 'FontSize', fs);
ylabel('b  (10^{-3} m^{-1})', 'FontSize', fs);
legend('Location','northeast','FontSize',fs,'Box','off');
title('log_{10}(J) — cost surface', 'FontWeight','normal','FontSize',fs);
set(gca,'FontSize',fs,'Box','off');

fig_dir = fullfile(fileparts(mfilename('fullpath')), '..', '..', 'docs', 'figures');
if ~exist(fig_dir,'dir'), mkdir(fig_dir); end
saveas(gcf, fullfile(fig_dir, 'toy_cost_surface.png'));

% ---------------------------------------------------------------
% 8. Figure 2: data vs fit
% ---------------------------------------------------------------

z_plot  = linspace(0, 550, 200);
P_prior_line = P0_prior * exp(-b_prior * z_plot);
P_fit_line   = P0_fit   * exp(-b_fit   * z_plot);
P_true_line  = P0_true  * exp(-b_true  * z_plot);

figure('Units','centimeters','Position',[13 2 8 9],'Color','white');
hold on;
errorbar(P_obs, z_obs, [], [], sigma_d, sigma_d, ...
    'bo', 'MarkerSize', 4, 'LineWidth', 0.8, 'DisplayName', 'obs');
plot(P_prior_line, z_plot, 'b--', 'LineWidth', 1.0, 'DisplayName', 'prior');
plot(P_true_line,  z_plot, 'k:',  'LineWidth', 1.2, 'DisplayName', 'true');
plot(P_fit_line,   z_plot, 'r-',  'LineWidth', 1.2, 'DisplayName', 'fit');
set(gca,'YDir','reverse','FontSize',fs,'Box','off');
xlabel('P (a.u.)', 'FontSize', fs);
ylabel('depth (m)', 'FontSize', fs);
legend('Location','southeast','FontSize',fs,'Box','off');
title('toy inverse: data vs fit', 'FontWeight','normal','FontSize',fs);

saveas(gcf, fullfile(fig_dir, 'toy_data_vs_fit.png'));
fprintf('Saved toy_cost_surface.png and toy_data_vs_fit.png\n');

% ---------------------------------------------------------------
% Helper functions
% ---------------------------------------------------------------

function P = fwd(p, z)
% forward model: exponential decay
P = p(1) * exp(-p(2) * z);
end

function J = cost_fn(p, z_obs, P_obs, sigma_d, P0_pr, b_pr, sig_P0, sig_b)
% cost = data misfit + prior penalty
P_mod   = fwd(p, z_obs);
misfit  = (P_obs - P_mod) ./ sigma_d;
J = sum(misfit.^2) + ((p(1)-P0_pr)/sig_P0)^2 + ((p(2)-b_pr)/sig_b)^2;
end
