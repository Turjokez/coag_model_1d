function sim = solve_column_1d_core(cfg, do_diff, do_coag, do_frag)
% solve_column_1d_core
% One shared 1-D column solver used by step scripts.

% --- grid and time
z_m = (0:cfg.dz_m:cfg.z_max_m)';
nz = numel(z_m);
t_s = (0:cfg.dt_s:cfg.t_max_s)';
nt = numel(t_s);

% --- sizes and speeds
size_um = cfg.size_um(:);
ns = numel(size_um);
speed_m_day = cfg.speed_m_s(:) .* 86400.0;

if isfield(cfg, 'speed_profile_m_s') && ~isempty(cfg.speed_profile_m_s)
    sp = cfg.speed_profile_m_s;
    if size(sp, 1) ~= nz || size(sp, 2) ~= ns
        error('speed_profile_m_s must be nz x ns');
    end
    speed_profile_m_s = sp;
else
    speed_profile_m_s = repmat(cfg.speed_m_s(:)', nz, 1);
end

if isfield(cfg, 'kz_profile_m2_s') && ~isempty(cfg.kz_profile_m2_s)
    kz_profile_m2_s = cfg.kz_profile_m2_s(:);
    if numel(kz_profile_m2_s) ~= nz
        error('kz_profile_m2_s must have nz points');
    end
elseif isfield(cfg, 'kz_m2_s') && ~isempty(cfg.kz_m2_s)
    kz_profile_m2_s = cfg.kz_m2_s .* ones(nz, 1);
else
    kz_profile_m2_s = zeros(nz, 1);
end

% --- initial condition
pulse_amp = cfg.pulse_amp;
if isscalar(pulse_amp)
    pulse_amp = pulse_amp .* ones(ns, 1);
else
    pulse_amp = pulse_amp(:);
end

top_mask = z_m <= 50.0;
conc = zeros(nz, ns, nt);
conc(top_mask, :, 1) = repmat(reshape(pulse_amp, 1, []), sum(top_mask), 1);

% --- per-particle volume
d_m = size_um .* 1e-6;
vol_part_m3 = (pi/6) .* (d_m .^ 3);

% --- tracked variables
column_number = zeros(nt, ns);
column_volume_by_size = zeros(nt, ns);
column_volume_total = zeros(nt, 1);
export_volume_total = zeros(nt, 1);
tracked_volume_total = zeros(nt, 1);
total_number = zeros(nt, 1);
export_number = zeros(ns, 1);
bottom_signal = zeros(nt, ns);

for is = 1:ns
    col_num = sum(conc(:, is, 1)) .* cfg.dz_m;
    column_number(1, is) = col_num;
    column_volume_by_size(1, is) = col_num .* vol_part_m3(is);
end
column_volume_total(1) = sum(column_volume_by_size(1, :));
tracked_volume_total(1) = column_volume_total(1);
total_number(1) = sum(column_number(1, :));
bottom_signal(1, :) = conc(end, :, 1);

% --- process settings
scheme = "upwind";
if isfield(cfg, 'scheme') && ~isempty(cfg.scheme)
    scheme = lower(string(cfg.scheme));
end

if isfield(cfg, 'coag_substeps') && ~isempty(cfg.coag_substeps)
    coag_substeps = max(1, round(cfg.coag_substeps));
else
    coag_substeps = 1;
end
dt_sub = cfg.dt_s ./ coag_substeps;

if isfield(cfg, 'frag_substeps') && ~isempty(cfg.frag_substeps)
    frag_substeps = max(1, round(cfg.frag_substeps));
else
    frag_substeps = 1;
end
dt_frag = cfg.dt_s ./ frag_substeps;

beta_m3_s = zeros(ns, ns);
if do_coag || do_frag
    beta_m3_s = local_build_beta_matrix(size_um .* 1e-4, cfg);
end

% for step-3 log
diff_alpha = max(kz_profile_m2_s) .* cfg.dt_s ./ max(cfg.dz_m .* cfg.dz_m, realmin);

% --- main time loop
for it = 2:nt
    c_prev = conc(:, :, it - 1);
    c_next = zeros(nz, ns);

    for is = 1:ns
        c_col = c_prev(:, is);
        w_col = speed_profile_m_s(:, is);

        % advection
        if scheme == "lax_wendroff"
            c_adv = local_lax_wendroff_step(c_col, w_col, cfg.dt_s, cfg.dz_m, 0.0);
        else
            c_adv = local_upwind_step(c_col, w_col, cfg.dt_s, cfg.dz_m, 0.0);
        end

        % diffusion
        if do_diff && any(kz_profile_m2_s > 0)
            c_new = local_diffusion_flux_step(c_adv, kz_profile_m2_s, cfg.dt_s, cfg.dz_m);
        else
            c_new = c_adv;
        end

        if scheme ~= "lax_wendroff"
            c_new(c_new < 0) = 0;
        end
        c_next(:, is) = c_new;

        % bottom export from previous state
        export_number(is) = export_number(is) + max(w_col(end), 0.0) .* max(c_col(end), 0.0) .* cfg.dt_s;
    end

    % coagulation
    if do_coag
        for iz = 1:nz
            n_vec = c_next(iz, :)';
            for sub = 1:coag_substeps
                dn = zeros(ns, 1);
                for i = 1:ns
                    for j = i:ns
                        if i == j
                            coll = beta_m3_s(i, j) .* n_vec(i) .* n_vec(j) .* dt_sub;
                            coll = min(coll, 0.25 .* n_vec(i));
                            k = min(ns, i + 1);
                            dn(i) = dn(i) - 2.0 .* coll;
                            dn(k) = dn(k) + coll;
                        else
                            coll = beta_m3_s(i, j) .* n_vec(i) .* n_vec(j) .* dt_sub;
                            coll = min(coll, 0.25 .* min(n_vec(i), n_vec(j)));
                            k = min(ns, max(i, j) + 1);
                            dn(i) = dn(i) - coll;
                            dn(j) = dn(j) - coll;
                            dn(k) = dn(k) + coll;
                        end
                    end
                end
                n_vec = n_vec + dn;
                n_vec(n_vec < 0) = 0;
            end
            c_next(iz, :) = n_vec';
        end
    end

    % fragmentation (simple conservative breakup)
    if do_frag
        c3 = 0.0;
        c4 = 1.0;
        if isfield(cfg, 'c3') && ~isempty(cfg.c3), c3 = cfg.c3; end
        if isfield(cfg, 'c4') && ~isempty(cfg.c4), c4 = cfg.c4; end
        eps_ref = 1.0;
        if isfield(cfg, 'epsilon_mks') && ~isempty(cfg.epsilon_mks)
            eps_ref = max(cfg.epsilon_mks, realmin);
        end

        frag_rate = c3 .* (size_um ./ 1000.0) .^ c4 .* sqrt(eps_ref);
        frag_rate(~isfinite(frag_rate)) = 0;
        frag_rate(frag_rate < 0) = 0;

        for iz = 1:nz
            n_vec = c_next(iz, :)';
            for sub = 1:frag_substeps
                dn = zeros(ns, 1);
                for i = 2:ns
                    loss = frag_rate(i) .* n_vec(i) .* dt_frag;
                    loss = min(loss, 0.30 .* n_vec(i));
                    dn(i) = dn(i) - loss;
                    gain = loss .* (vol_part_m3(i) ./ max(vol_part_m3(i - 1), realmin));
                    dn(i - 1) = dn(i - 1) + gain;
                end
                n_vec = n_vec + dn;
                n_vec(n_vec < 0) = 0;
            end
            c_next(iz, :) = n_vec';
        end
    end

    % store
    conc(:, :, it) = c_next;
    bottom_signal(it, :) = c_next(end, :);

    for is = 1:ns
        col_num = sum(c_next(:, is)) .* cfg.dz_m;
        column_number(it, is) = col_num;
        column_volume_by_size(it, is) = col_num .* vol_part_m3(is);
    end
    export_volume_total(it) = sum(export_number .* vol_part_m3);
    column_volume_total(it) = sum(column_volume_by_size(it, :));
    tracked_volume_total(it) = column_volume_total(it) + export_volume_total(it);
    total_number(it) = sum(column_number(it, :));
end

% --- output struct
sim = struct();
sim.cfg = cfg;
sim.t_s = t_s;
sim.z_m = z_m;
sim.size_um = size_um;
sim.speed_m_day = speed_m_day;
sim.speed_profile_m_s = speed_profile_m_s;
sim.speed_profile_m_day = speed_profile_m_s .* 86400.0;
sim.conc = conc;
sim.bottom_signal = bottom_signal;
sim.column_number = column_number;
sim.total_number = total_number;
sim.column_volume_by_size = column_volume_by_size;
sim.column_volume_total = column_volume_total;
sim.export_volume_total = export_volume_total;
sim.tracked_volume_total = tracked_volume_total;
sim.tracked_mass = tracked_volume_total; % keep old script compatibility
sim.beta_m3_s = beta_m3_s;
sim.diff_alpha = diff_alpha;
sim.cfl = struct('max_cfl', max(abs(speed_profile_m_s(:))) .* cfg.dt_s ./ cfg.dz_m);
end

function c_next = local_upwind_step(c_prev, w_col, dt_s, dz_m, c_in)
flux = local_adv_flux_upwind(c_prev, w_col, c_in);
c_next = c_prev - (dt_s ./ dz_m) .* (flux(2:end) - flux(1:end-1));
end

function c_next = local_lax_wendroff_step(c_prev, w_col, dt_s, dz_m, c_in)
% Simple scalar LW using local mean speed. Used only for compare checks.
n = numel(c_prev);
c = c_prev(:);
if isscalar(w_col)
    u = w_col .* ones(n, 1);
else
    u = w_col(:);
end

% inlet
c_ext = [c_in; c; c(end)];
u_ext = [u(1); u; u(end)];
c_next = c;
for j = 1:n
    uj = u_ext(j + 1);
    r = uj .* dt_s ./ dz_m;
    cjm = c_ext(j);
    cj = c_ext(j + 1);
    cjp = c_ext(j + 2);
    c_next(j) = cj - 0.5 .* r .* (cjp - cjm) + 0.5 .* r .* r .* (cjp - 2 .* cj + cjm);
end
end

function c_next = local_diffusion_flux_step(c_prev, kz_m2_s, dt_s, dz_m)
flux = local_diff_flux(c_prev, kz_m2_s, dz_m);
c_next = c_prev - (dt_s ./ dz_m) .* (flux(2:end) - flux(1:end-1));
end

function flux = local_adv_flux_upwind(c_prev, w_m_s, c_in)
c_prev = c_prev(:);
n = numel(c_prev);
if isscalar(w_m_s)
    w_col = w_m_s .* ones(n, 1);
else
    w_col = w_m_s(:);
end
flux = zeros(n + 1, 1);
flux(1) = max(w_col(1), 0.0) .* c_in;
for j = 1:(n - 1)
    w_face = 0.5 .* (w_col(j) + w_col(j + 1));
    if w_face >= 0
        c_up = c_prev(j);
    else
        c_up = c_prev(j + 1);
    end
    flux(j + 1) = w_face .* c_up;
end
flux(n + 1) = max(w_col(end), 0.0) .* c_prev(end);
end

function flux = local_diff_flux(c_prev, kz_m2_s, dz_m)
c_prev = c_prev(:);
kz_m2_s = kz_m2_s(:);
n = numel(c_prev);
flux = zeros(n + 1, 1);
flux(1) = 0.0;
for j = 1:(n - 1)
    k_face = 0.5 .* (kz_m2_s(j) + kz_m2_s(j + 1));
    grad = (c_prev(j + 1) - c_prev(j)) ./ dz_m;
    flux(j + 1) = -k_face .* grad;
end
flux(n + 1) = 0.0;
end

function beta_m3_s = local_build_beta_matrix(size_cm, cfg)
[D1, D2] = ndgrid(size_cm, size_cm);

mode_name = "both";
if isfield(cfg, 'kernel_mode') && ~isempty(cfg.kernel_mode)
    mode_name = lower(string(cfg.kernel_mode));
end

[beta_ds_cm3_s, ~, ~] = local_beta_diff_sed_from_law(D1, D2, cfg.law_name);
beta_ds_m3_s = beta_ds_cm3_s .* 1e-6;

eps_mks = 1e-6;
if isfield(cfg, 'epsilon_mks') && ~isempty(cfg.epsilon_mks)
    eps_mks = cfg.epsilon_mks;
end
rg_m = 0.5 .* (D1 + D2) .* 1e-2;
beta_shear_m3_s = sqrt(max(eps_mks, 0)) .* (rg_m .^ 3);

switch mode_name
    case "shear_only"
        beta_m3_s = beta_shear_m3_s;
    case "diff_sed_only"
        beta_m3_s = beta_ds_m3_s;
    otherwise
        beta_m3_s = beta_shear_m3_s + beta_ds_m3_s;
end

if isfield(cfg, 'scale_shear') && ~isempty(cfg.scale_shear) && mode_name ~= "diff_sed_only"
    beta_m3_s = beta_m3_s .* cfg.scale_shear;
end
if isfield(cfg, 'scale_diff_sed') && ~isempty(cfg.scale_diff_sed) && mode_name ~= "shear_only"
    beta_m3_s = beta_m3_s .* cfg.scale_diff_sed;
end
if isfield(cfg, 'coag_scale') && ~isempty(cfg.coag_scale)
    beta_m3_s = beta_m3_s .* cfg.coag_scale;
end

beta_m3_s(~isfinite(beta_m3_s)) = 0;
beta_m3_s(beta_m3_s < 0) = 0;
end

function [beta, w1, w2] = local_beta_diff_sed_from_law(d1_cm, d2_cm, law_name)
w1 = sinking_speed_named(d1_cm, law_name);
w2 = sinking_speed_named(d2_cm, law_name);
beta = (pi/4) .* (d1_cm + d2_cm) .* (d1_cm + d2_cm) .* abs(w1 - w2);
end

