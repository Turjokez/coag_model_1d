Directory and File Structure
============================

The directory tree of ``1d-model-testing/`` is::

    1d-model-testing/
      src/                       % MATLAB class definitions
        SimulationConfig.m       % all tunable parameters
        ColumnGrid.m             % depth grid (n_z layers, dz spacing)
        ColumnSimulation.m       % top-level run object
        ColumnRHS.m              % right-hand side: all physics
        ColumnTransport.m        % first-order upwind transport
        ZooplanktonGrazing.m     % Stemmann 2004 grazing model
        Disaggregation.m         % operator-split disaggregation
        FecalCrossCoag.m         % fecal--aggregate cross-coagulation
        DepthProfile.m           % depth profiles of eps, T, S, zoo
        KernelLibrary.m          % Brownian, shear, DS collision kernels
        SettlingVelocityService.m  % sinking speed laws

      scripts/
        data/                    % data loaders and main run scripts
        examples/                % three worked examples (Section 6)
          run_example_01.m
          run_example_02.m
          run_example_03.m
        inverse_toy/             % toy inverse problem (Section 8)
          run_toy_inverse.m

      data/
        NA/
          Turbulance/keps_for_dave.mat   % daily turbulence profiles
          uvp/raw/                       % UVP .sb observation files

      docs/
        user_manual/             % LaTeX version of this document
        sphinx_manual/           % Sphinx version of this document
        figures/                 % figures saved by example scripts

The class ``SimulationConfig`` is the central configuration object.
All parameters that control the model physics, numerics, and biological
processes are properties of this class with documented default values.
The complete parameter list is in :doc:`parameters`.
