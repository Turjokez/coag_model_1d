function J = cost_fn_col(params, obs, obs_depths, col_grid, keps_day, prof, phi_bc_daily, n_days, cfg_base, prior)
% cost_fn_col  Cost function for column model inverse fitting.
%
% params = [alpha, zoo_c_scale]
%
% Cost = data misfit (log-space) + prior penalty
%
% Log-space misfit is used because BV spans orders of magnitude
% and observation uncertainty is roughly proportional to the value.

alpha_try = params(1);
zoo_c_sc  = params(2);

% --- bounds check: return large J if params out of range ---
if alpha_try <= 0 || alpha_try > 2 || zoo_c_sc <= 0 || zoo_c_sc > 10
    J = 1e6;
    return
end

% --- run forward model ---
phi_mod = fwd_column(params, obs_depths, col_grid, keps_day, prof, phi_bc_daily, n_days, cfg_base);
bv_mod  = sum(phi_mod, 2);   % total BV at each obs depth [n_dep x 1]

% --- log-space misfit ---
% sigma_log = 0.5 in natural log (~ factor-of-1.6 obs uncertainty)
sigma_log = 0.5;
J_data = 0;
for id = 1:numel(obs_depths)
    if bv_mod(id) <= 0
        J = 1e6; return
    end
    resid   = log(bv_mod(id)) - log(obs.bv_total(id));
    J_data  = J_data + (resid / sigma_log)^2;
end

% --- prior penalty ---
% alpha: log-normal prior, mean=prior.alpha, sigma=prior.sigma_alpha
resid_alpha = (log(alpha_try) - log(prior.alpha)) / prior.sigma_log_alpha;

% zoo_c_scale: log-normal prior, mean=1 (no scaling), sigma=prior.sigma_log_zoo
resid_zoo   = (log(zoo_c_sc) - log(prior.zoo_c_scale)) / prior.sigma_log_zoo;

J_prior = resid_alpha^2 + resid_zoo^2;

J = J_data + J_prior;

fprintf('  alpha=%.3f  zoo_c_sc=%.2f  J_data=%.2f  J_prior=%.2f  J=%.2f\n', ...
    alpha_try, zoo_c_sc, J_data, J_prior, J);
end
