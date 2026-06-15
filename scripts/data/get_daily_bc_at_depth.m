function bc = get_daily_bc_at_depth(uvp_file, cfg, col_grid, bc_depth_m, k_compare)
% GET_DAILY_BC_AT_DEPTH  Build daily UVP boundary condition at one depth.
%
% Usage:
%   bc = get_daily_bc_at_depth(uvp_file, cfg, col_grid, 100, 3:10)
%
% Output:
%   bc.dates          - daily dates [YYYYMMDD]
%   bc.n_days         - number of days
%   bc.phi_bc_daily   - daily BC in model bins [n_days x n_sec]
%   bc.uvpd           - parsed UVP daily struct
%   bc.iz_bc          - UVP depth index used for BC
%   bc.bc_depth_m     - actual UVP depth used
%   bc.d_model_um     - model bin diameters [um]
%   bc.mask_uvp_model - model bins in 100-2000 um
%   bc.mask_small     - model bins below 100 um
%   bc.id_model_best  - best day index in model daily series
%   bc.id_uvp_best    - best day index in UVP daily series
%   bc.best_date      - best cast date [YYYYMMDD]
%   bc.phi_uvp_cmp    - UVP phi at comparison depths [n_cmp x 1]
%   bc.phi_uvp_spec   - UVP phi by bin at comparison depths [n_cmp x n_bins]
%   bc.d_uvp_ok       - UVP bin centers in 100-2000 um
%   bc.dw_uvp_ok      - UVP bin widths in 100-2000 um
%
% Notes:
%   - UVP range 100-2000 um is mapped by overlap fraction.
%   - Model bins below 100 um are filled by a log-log power law fit
%     to the first UVP bins from 100 to 400 um.

uvpd = parse_uvp_daily(uvp_file);
[~, iz_bc] = min(abs(uvpd.depth_m - bc_depth_m));

grid       = cfg.derive();
n_sec      = cfg.n_sections;
r_cm       = (0.75 / pi * grid.av_vol(:)).^(1/3);
d_model_um = 2 * r_cm * 1e4;

mask_uvp_data  = uvpd.d_um >= 100 & uvpd.d_um < 2000;
mask_small     = d_model_um < 100;
mask_uvp_model = d_model_um >= 100 & d_model_um < 2000;

d_uvp_ok  = uvpd.d_um(mask_uvp_data);
dw_uvp_ok = uvpd.dw(mask_uvp_data);
n_uvp_ok  = sum(mask_uvp_data);

% model bin edges from geometric midpoint
d_model_edges = zeros(1, n_sec + 1);
d_model_edges(1)       = d_model_um(1)^2 / d_model_um(2);
d_model_edges(n_sec+1) = d_model_um(n_sec)^2 / d_model_um(n_sec-1);
for k = 2:n_sec
    d_model_edges(k) = sqrt(d_model_um(k-1) * d_model_um(k));
end

% overlap mapping for UVP-visible bins
overlap_frac = zeros(n_sec, n_uvp_ok);
for j = 1:n_uvp_ok
    uvp_lo = d_uvp_ok(j) - dw_uvp_ok(j);
    uvp_hi = d_uvp_ok(j);
    for k = 1:n_sec
        lo = max(d_model_edges(k),   uvp_lo);
        hi = min(d_model_edges(k+1), uvp_hi);
        if hi > lo
            overlap_frac(k, j) = (hi - lo) / dw_uvp_ok(j);
        end
    end
end

% use full date span, fill gaps by nearest day if needed
dn_obs   = datenum(num2str(uvpd.dates), 'yyyymmdd');
dn_all   = (dn_obs(1):dn_obs(end))';
dates_all = str2double(cellstr(datestr(dn_all, 'yyyymmdd')));
n_days    = numel(dn_all);

phi_bc_daily = zeros(n_days, n_sec);

for id_m = 1:n_days
    [is_hit, ib] = ismember(dates_all(id_m), uvpd.dates);
    if is_hit
        id_u = ib;
    else
        [~, id_u] = min(abs(dn_obs - dn_all(id_m)));
    end

    phi_row = squeeze(uvpd.phi(id_u, iz_bc, mask_uvp_data));
    phi_row(isnan(phi_row)) = 0;

    phi_mod = zeros(1, n_sec);

    % visible UVP bins -> model bins
    for k = 1:n_sec
        phi_mod(k) = sum(overlap_frac(k, :) .* phi_row(:)');
    end

    % power-law fill below 100 um
    fit_ok = d_uvp_ok >= 100 & d_uvp_ok <= 400 & phi_row(:)' > 0;
    if sum(fit_ok) >= 2
        p = polyfit(log10(d_uvp_ok(fit_ok)), log10(phi_row(fit_ok)'), 1);
        phi_fill = 10 .^ polyval(p, log10(d_model_um(mask_small)));
        phi_fill(~isfinite(phi_fill) | phi_fill < 0) = 0;
        phi_mod(mask_small) = phi_fill;
    end

    phi_bc_daily(id_m, :) = phi_mod;
end

% best day = strongest BC in UVP-visible range
phi_vis = sum(phi_bc_daily(:, mask_uvp_model), 2);
[~, id_model_best] = max(phi_vis);
[~, id_uvp_best] = min(abs(dn_obs - dn_all(id_model_best)));

z_compare = col_grid.z_centers(k_compare);
phi_uvp_cmp = zeros(numel(k_compare), 1);
phi_uvp_spec = zeros(numel(k_compare), n_uvp_ok);
for i = 1:numel(k_compare)
    [~, iz_u] = min(abs(uvpd.depth_m - z_compare(i)));
    phi_row = squeeze(uvpd.phi(id_uvp_best, iz_u, mask_uvp_data));
    phi_row(isnan(phi_row)) = 0;
    phi_uvp_cmp(i) = sum(phi_row);
    phi_uvp_spec(i, :) = phi_row;
end

bc.dates          = dates_all;
bc.n_days         = n_days;
bc.phi_bc_daily   = phi_bc_daily;
bc.uvpd           = uvpd;
bc.iz_bc          = iz_bc;
bc.bc_depth_m     = uvpd.depth_m(iz_bc);
bc.d_model_um     = d_model_um;
bc.mask_uvp_model = mask_uvp_model;
bc.mask_small     = mask_small;
bc.id_model_best  = id_model_best;
bc.id_uvp_best    = id_uvp_best;
bc.best_date      = dates_all(id_model_best);
bc.phi_uvp_cmp    = phi_uvp_cmp;
bc.phi_uvp_spec   = phi_uvp_spec;
bc.z_compare      = z_compare;
bc.d_uvp_ok       = d_uvp_ok;
bc.dw_uvp_ok      = dw_uvp_ok;

end
