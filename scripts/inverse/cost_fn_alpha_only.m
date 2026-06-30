function J = cost_fn_alpha_only(alpha_try, obs, obs_depths, col_grid, keps_day, prof, phi_bc_daily, n_days, cfg_base, prior)
% cost_fn_alpha_only  Cost for alpha-only fit at deep depths.
%
% alpha_try = scalar alpha
% zoo_c stays fixed at the base Stemmann value.

if alpha_try <= 0 || alpha_try > 2
    J = 1e6;
    return
end

phi_mod = fwd_column([alpha_try, 1.0], obs_depths, col_grid, keps_day, prof, phi_bc_daily, n_days, cfg_base);
bv_mod  = sum(phi_mod, 2);

sigma_log = 0.5;
J_data = 0;
for id = 1:numel(obs_depths)
    if bv_mod(id) <= 0
        J = 1e6;
        return
    end
    resid  = log(bv_mod(id)) - log(obs.bv_total(id));
    J_data = J_data + (resid / sigma_log)^2;
end

resid_alpha = (log(alpha_try) - log(prior.alpha)) / prior.sigma_log_alpha;
J_prior = resid_alpha^2;

J = J_data + J_prior;

fprintf('  alpha=%.4f  J_data=%.2f  J_prior=%.2f  J=%.2f\n', ...
    alpha_try, J_data, J_prior, J);
end
