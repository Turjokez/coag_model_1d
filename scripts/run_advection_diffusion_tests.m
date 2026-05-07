% run_advection_diffusion_tests
% Steps:
% 1. run advection-diffusion cases
% 2. validate spread and stability
% 3. save short logs

clear;
clc;

repo_root = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(fullfile(repo_root, 'src')));

fig_dir = fullfile(repo_root, 'output', 'figures');
log_dir = fullfile(repo_root, 'output', 'logs');
table_dir = fullfile(repo_root, 'output', 'tables');

if ~exist(fig_dir, 'dir')
    mkdir(fig_dir);
end
if ~exist(log_dir, 'dir')
    mkdir(log_dir);
end
if ~exist(table_dir, 'dir')
    mkdir(table_dir);
end

% Keep the same baseline settings so the new spread is easier to read.
law_name = 'kriest_8';
scheme = 'lax_wendroff';
size_um = [100; 500; 1000; 3000];
size_cm = size_um ./ 1e4;
speed_cm_s = local_named_speed(size_cm, law_name);
speed_m_s = speed_cm_s .* 0.01;

z_max_m = 1000.0;
dz_m = 5.0;
dt_s = 0.90 .* dz_m ./ max(speed_m_s);
travel_s = z_max_m ./ speed_m_s;
t_max_s = 1.20 .* max(travel_s);
kz_m2_s = 1e-4;

cfg = struct();
cfg.z_max_m = z_max_m;
cfg.dz_m = dz_m;
cfg.dt_s = dt_s;
cfg.t_max_s = t_max_s;
cfg.size_um = size_um;
cfg.speed_m_s = speed_m_s;
cfg.pulse_amp = 1.0;
cfg.law_name = law_name;
cfg.scheme = scheme;

sim_adv = solve_advection_only(cfg);

cfg_diff = cfg;
cfg_diff.kz_m2_s = kz_m2_s;
sim_diff = solve_advection_diffusion(cfg_diff);

val = validate_diffusion(sim_adv, sim_diff, fig_dir, log_dir, table_dir);

disp('Saved diffusion log:');
disp(val.log_path);
disp('Saved diffusion table:');
disp(val.csv_path);
disp('Summary:');
disp(val.summary);

function w = local_named_speed(diam_cm, law_name)
switch lower(string(law_name))
    case "current"
        w = sinking_speed_current(diam_cm);
    case "kriest_8"
        w = sinking_speed_kriest8(diam_cm);
    case "kriest_9"
        w = sinking_speed_kriest9(diam_cm);
    case "siegel_2025"
        w = sinking_speed_siegel2025(diam_cm);
    otherwise
        error('run_advection_diffusion_tests:law', 'Unknown law: %s', law_name);
end
end
