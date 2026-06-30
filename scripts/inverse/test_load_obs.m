% test_load_obs.m
% Load UVP observations at 125, 325, 475 m and print them.
% Check the depth levels used and BV values before building cost function.

script_dir = fileparts(mfilename('fullpath'));
addpath(script_dir);
addpath(fullfile(script_dir, '..', '..', 'src'));
addpath(fullfile(script_dir, '..', 'data'));

uvp_file = fullfile(script_dir, '..', '..', 'data', 'NA', 'uvp', 'raw', ...
    '745f00117e_EXPORTS-EXPORTSNA_UVP5-ParticulateLevel2_differential_survey_20210504-20210529_R1.sb');

obs_depths = [125, 325, 475];
obs = load_uvp_obs(uvp_file, obs_depths);

fprintf('UVP observations (100-2000 um, averaged over all casts):\n');
for id = 1:numel(obs_depths)
    fprintf('  target %d m -> actual %d m: BV = %.4e cm^3/cm^3\n', ...
        obs_depths(id), obs.depth(id), obs.bv_total(id));
end
