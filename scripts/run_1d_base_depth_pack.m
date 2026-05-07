% run_1d_base_depth_pack
% Short note:
% 1. run base 1-D depth checks in order
% 2. check key output files exist
% 3. save one simple pass/fail log

clear;
clc;

repo_root = fileparts(fileparts(mfilename('fullpath')));
script_dir = fullfile(repo_root, 'scripts');
out_fig = fullfile(repo_root, 'output', 'figures');
out_tab = fullfile(repo_root, 'output', 'tables');
out_log = fullfile(repo_root, 'output', 'logs');
matlab_bin = '/Applications/MATLAB_R2025a.app/bin/matlab';

if ~exist(out_log, 'dir')
    mkdir(out_log);
end

steps = {
    'run_1d_base_depth_step1.m'
};

step_ok = false(numel(steps), 1);
step_msg = strings(numel(steps), 1);

for i = 1:numel(steps)
    this_file = fullfile(script_dir, steps{i});
    try
        cmd = sprintf('"%s" -batch "cd(''%s''); run(''scripts/%s'');"', ...
            matlab_bin, repo_root, steps{i});
        [status, out_txt] = system(cmd);
        if status == 0
            step_ok(i) = true;
            step_msg(i) = "ok";
        else
            step_ok(i) = false;
            step_msg(i) = "failed: " + string(strtrim(out_txt));
        end
    catch ME
        step_ok(i) = false;
        step_msg(i) = "failed: " + string(ME.message);
    end
end

need_files = {
    fullfile(out_fig, 'base_depth_step1_pulse_profiles.png')
    fullfile(out_fig, 'base_depth_step1_depth_size_snapshots.png')
    fullfile(out_fig, 'base_depth_step1_conservation.png')
    fullfile(out_tab, 'base_depth_step1_speed_summary.csv')
    fullfile(out_log, 'base_depth_step1_check.txt')
};

file_ok = false(numel(need_files), 1);
for i = 1:numel(need_files)
    file_ok(i) = isfile(need_files{i});
end

pack_log = fullfile(out_log, 'step_base_depth_pack.txt');
fid = fopen(pack_log, 'w');
fprintf(fid, 'Base 1-D depth pack check\n\n');
fprintf(fid, 'Run date: %s\n\n', datestr(now));

fprintf(fid, 'Step run status:\n');
for i = 1:numel(steps)
    fprintf(fid, '- %s : %s\n', steps{i}, step_msg(i));
end
fprintf(fid, '\n');

fprintf(fid, 'Required output files:\n');
for i = 1:numel(need_files)
    if file_ok(i)
        txt = 'ok';
    else
        txt = 'missing';
    end
    fprintf(fid, '- %s : %s\n', need_files{i}, txt);
end
fprintf(fid, '\n');

if all(step_ok) && all(file_ok)
    fprintf(fid, 'Final status: PASS\n');
else
    fprintf(fid, 'Final status: NEEDS WORK\n');
end
fclose(fid);

disp('Saved pack log:');
disp(pack_log);
