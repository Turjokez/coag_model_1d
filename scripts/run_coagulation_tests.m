% run_coagulation_tests
% Short note:
% 1. keep one entry point for the coagulation stage
% 2. use the trusted step-4 setup

this_dir = fileparts(mfilename('fullpath'));
run(fullfile(this_dir, 'run_1d_step4_coagulation_check.m'));
