Obtaining and Installing the Code
==================================

Requirements
------------

The model requires MATLAB version R2021a or later. No additional toolboxes
beyond the core MATLAB installation are required. The code has been tested
on macOS and Linux.

Directory Layout
----------------

The code is organized into two directories that must reside at the same
level on the user's computer::

    coag_model/
        coag_model_final/    % shared coagulation kernel library
        1d-model-testing/    % 1-D column model scripts, data, and documentation

The directory ``coag_model_final/`` contains the sectional coagulation solver
shared between the 1-D model and earlier slab-model tests.
The directory ``1d-model-testing/`` contains all scripts, source classes,
observational data, and documentation specific to the 1-D column model.

Setting MATLAB Paths
--------------------

Each script sets its own MATLAB path. The pattern at the top of every script is:

.. code-block:: matlab

    script_dir = fileparts(mfilename('fullpath'));
    addpath(script_dir);
    addpath(fullfile(script_dir, '..', '..', 'src'));

The first ``addpath`` call adds the script's own directory so helper functions
in the same folder are found. The second adds the ``src/`` directory, which
contains all model class definitions. Do not add paths globally across MATLAB
sessions; the per-script pattern avoids conflicts between experiments.
