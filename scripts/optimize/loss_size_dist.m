function L = loss_size_dist(Y_model, Y_obs)
% LOSS_SIZE_DIST  Sum of squared differences in size spectrum.
%
% Y_model : (n_z x n_sec) biovolume array from model
% Y_obs   : (n_z x n_sec) observed (or target) biovolume array
% L       : scalar loss (sum of squared differences)
%
% Note: log-space comparison weights all bins equally.
% Linear-space would make large bins dominate.
% Replace Y_obs with EXPORTS UVP data when available.

eps_floor = 1e-30;   % floor to avoid log(0)

log_model = log10(Y_model + eps_floor);
log_obs   = log10(Y_obs   + eps_floor);

L = sum((log_model - log_obs).^2, 'all');
end
