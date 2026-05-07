% run_fragmentation_tests
% Short note:
% 1. keep one entry point for the fragmentation stage
% 2. use the trusted step-5 setup

this_dir = fileparts(mfilename('fullpath'));
run(fullfile(this_dir, 'run_1d_step5_fragmentation_check.m'));
