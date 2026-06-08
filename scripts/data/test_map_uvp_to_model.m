% test_map_uvp_to_model.m
% Check that UVP data can become a model Y0.

script_dir = fileparts(mfilename('fullpath'));
addpath(script_dir);
addpath(fullfile(script_dir, '..', '..', 'src'));

sb_file = fullfile(script_dir, '..', '..', 'data', 'NA', 'uvp', 'raw', ...
    '745f00117e_EXPORTS-EXPORTSNA_UVP5-ParticulateLevel2_differential_survey_20210504-20210529_R1.sb');

fprintf('=== test_map_uvp_to_model ===\n');

uvp = parse_uvp(sb_file);

cfg = SimulationConfig('n_sections', 30);
col_grid = ColumnGrid(1000, 20);

mapped = map_uvp_to_model(uvp, cfg, col_grid);

pass = true;

% 1. Y0 shape
if all(size(mapped.Y0_surface) == [col_grid.n_z, cfg.n_sections])
    fprintf('PASS  Y0 size is %d x %d\n', col_grid.n_z, cfg.n_sections);
else
    fprintf('FAIL  wrong Y0 size\n');
    pass = false;
end

% 2. surface has mass
surface_sum = sum(mapped.Y0_surface(1, :));
if surface_sum > 0
    fprintf('PASS  surface phi sum = %.2e cm^3/cm^3\n', surface_sum);
else
    fprintf('FAIL  surface phi is zero\n');
    pass = false;
end

% 3. deeper layers are zero in Y0
deep_sum = sum(mapped.Y0_surface(2:end, :));
deep_sum = sum(deep_sum);
if deep_sum == 0
    fprintf('PASS  only surface layer is initialized\n');
else
    fprintf('FAIL  deep layers are not zero\n');
    pass = false;
end

% 4. UVP bins mapped inside model bins
if all(mapped.bin_map >= 1) && all(mapped.bin_map <= cfg.n_sections)
    fprintf('PASS  UVP bins mapped to model bins %d-%d\n', ...
        min(mapped.bin_map), max(mapped.bin_map));
else
    fprintf('FAIL  bin map out of range\n');
    pass = false;
end

if pass
    fprintf('ALL PASS\n');
else
    fprintf('SOME CHECKS FAILED\n');
end

% simple plot
figure;
bar(1:cfg.n_sections, mapped.surface_phi);
xlabel('section');
ylabel('\phi  [cm^3 cm^{-3}]');
title('UVP surface forcing mapped to model bins');
