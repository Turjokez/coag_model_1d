% run_1d_step5_fragmentation_check
% Steps:
% 1. keep the trusted upwind + diffusion + coag baseline
% 2. add fragmentation only
% 3. save simple figures and one short summary

clear;
clc;

repo_root = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(fullfile(repo_root, 'src')));

fig_dir = fullfile(repo_root, 'output', 'figures');
log_dir = fullfile(repo_root, 'output', 'logs');
tab_dir = fullfile(repo_root, 'output', 'tables');

if ~exist(fig_dir, 'dir')
    mkdir(fig_dir);
end
if ~exist(log_dir, 'dir')
    mkdir(log_dir);
end
if ~exist(tab_dir, 'dir')
    mkdir(tab_dir);
end

% Keep this first fragmentation check close to step 4.
law_name = 'kriest_8';
size_um = round(logspace(log10(200), log10(3000), 8))';
size_cm = size_um .* 1e-4;
pulse_amp = powerlaw_concentration(size_cm, 5e-3, -2.5);
speed_cm_s = sinking_speed_named(size_cm, law_name);
speed_m_s = speed_cm_s .* 0.01;

cfg = struct();
cfg.z_max_m = 1000.0;
cfg.dz_m = 5.0;
cfg.dt_s = 0.50 .* cfg.dz_m ./ max(speed_m_s);
cfg.t_max_s = 1.10 .* max(cfg.z_max_m ./ speed_m_s);
cfg.size_um = size_um;
cfg.speed_m_s = speed_m_s;
cfg.pulse_amp = pulse_amp;
cfg.law_name = law_name;
cfg.scheme = 'upwind';
cfg.kz_m2_s = 1e-4;
cfg.kernel_mode = 'shear_only';
cfg.epsilon_mks = 1e-6;
cfg.coag_scale = 100.0;
cfg.coag_substeps = 4;
cfg.scale_shear = 1.0;
cfg.scale_diff_sed = 0.0;

sim_coag = solve_with_coagulation(cfg);

cfg_frag = cfg;
cfg_frag.frag_substeps = 4;
cfg_frag.c3 = 0.005;
cfg_frag.c4 = 1.45;
sim_frag = solve_with_fragmentation(cfg_frag);

tag = 'step5_fragmentation';
val = validate_fragmentation_conservation(sim_coag, sim_frag, fig_dir, log_dir, tab_dir, tag);

disp('Saved step 5 figures and summary:');
disp(fullfile(fig_dir, [tag '_conservation.png']));
disp(fullfile(fig_dir, [tag '_small_size_volume.png']));
disp(fullfile(fig_dir, [tag '_column_psd.png']));
disp(fullfile(fig_dir, [tag '_final_depth_psd.png']));
disp(fullfile(fig_dir, [tag '_bottom_psd_time.png']));
disp(val.csv_path);
disp(val.log_path);
