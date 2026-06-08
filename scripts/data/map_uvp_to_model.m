function out = map_uvp_to_model(uvp, cfg, col_grid)
% MAP_UVP_TO_MODEL  Put UVP biovolume into model size bins.
%
% Usage:
%   out = map_uvp_to_model(uvp, cfg, col_grid)
%
% Output:
%   out.phi_depth  - UVP phi mapped to model bins, depth x section
%   out.Y0_surface - initial condition with top layer from UVP surface phi
%   out.d_model_um - model diameter [um]
%   out.bin_map    - UVP bin -> model bin index
%
% Notes:
%   - UVP phi is already [cm^3/cm^3].
%   - We add each UVP bin to the nearest model bin by diameter.
%   - Only the top 5 m mean is used for Y0_surface.

grid = cfg.derive();
n_sec = cfg.n_sections;
n_z = col_grid.n_z;

% model equivalent diameter from conserved volume
r_cm = (0.75 / pi * grid.av_vol(:)).^(1/3);
d_model_um = 2 * r_cm * 1e4;   % cm -> um

% nearest model bin for each UVP bin
bin_map = zeros(size(uvp.d_um));
for i = 1:numel(uvp.d_um)
    [~, bin_map(i)] = min(abs(d_model_um - uvp.d_um(i)));
end

% map full UVP depth profile
phi_depth = zeros(numel(uvp.depth_m), n_sec);
for i = 1:numel(uvp.d_um)
    k = bin_map(i);
    vals = uvp.phi(:, i);
    vals(isnan(vals)) = 0;
    phi_depth(:, k) = phi_depth(:, k) + vals;
end

% use top 5 m as model surface initial condition
surface_rows = uvp.depth_m <= 5;
if ~any(surface_rows)
    error('No UVP data found in top 5 m.');
end

surface_phi = mean_no_nan(phi_depth(surface_rows, :), 1);
surface_phi(isnan(surface_phi)) = 0;

Y0 = zeros(n_z, n_sec);
Y0(1, :) = surface_phi;

out.phi_depth  = phi_depth;
out.Y0_surface = Y0;
out.d_model_um = d_model_um(:)';
out.bin_map    = bin_map;
out.surface_phi = surface_phi;
end

function y = mean_no_nan(x, dim)
good = ~isnan(x);
x(~good) = 0;
n = sum(good, dim);
y = sum(x, dim) ./ max(n, 1);
y(n == 0) = NaN;
end
