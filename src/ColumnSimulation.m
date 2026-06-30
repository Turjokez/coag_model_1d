%> @brief Top-level driver for the 1-D depth column simulation.
%> @details Integrates particle biovolume forward in time over a vertical
%>          column using ColumnRHS for transport and process rates.
%> @par Example
%> @code
%>   cfg  = SimulationConfig('sinking_law','kriest_8','n_sections',30);
%>   grid = ColumnGrid(1000, 20);
%>   prof = DepthProfile.typical(grid.z_centers);
%>   sim  = ColumnSimulation(cfg, grid, prof);
%>   res  = sim.run();
%> @endcode
classdef ColumnSimulation < handle
    % COLUMNSIMULATION  Top-level 1-D depth column simulation.
    %
    % Analogous to CoagulationSimulation but for the full 1-D depth model.
    % Uses ColumnRHS for transport + process rates.
    %
    % Usage:
    %   cfg  = SimulationConfig();
    %   cfg.sinking_law  = 'kriest_8';
    %   cfg.n_sections   = 20;
    %   cfg.enable_coag  = true;
    %   cfg.disagg_mode  = 'operator_split';
    %   cfg.disagg_dmax_cm = 2.0;
    %
    %   cgrid  = ColumnGrid(2000, 40);           % 2000 m, 40 layers
    %   prof   = DepthProfile.typical(cgrid.z_centers);
    %   sim    = ColumnSimulation(cfg, cgrid, prof);
    %   result = sim.run();
    %
    % result fields:
    %   result.time          - t_output vector [day]
    %   result.concentrations - (n_t x n_z x n_sec) array
    %   result.col_grid      - ColumnGrid
    %   result.profile       - DepthProfile
    %   result.cfg           - SimulationConfig

    properties
        cfg        % SimulationConfig
        col_grid   % ColumnGrid
        profile    % DepthProfile
        size_grid  % DerivedGrid
        rhs        % ColumnRHS
    end

    methods
        %> @brief Constructor.
        %> @param cfg       SimulationConfig object.
        %> @param col_grid  ColumnGrid defining depth layers.
        %> @param profile   DepthProfile with epsilon, shear, zoo abundance vs depth.
        function obj = ColumnSimulation(cfg, col_grid, profile)
            obj.cfg       = cfg;
            obj.col_grid  = col_grid;
            obj.profile   = profile;
            obj.size_grid = DerivedGrid(cfg);
            obj.rhs       = ColumnRHS(cfg, obj.size_grid, col_grid, profile);
        end

        %> @brief Integrate the column forward in time.
        %> @param varargin  Optional: 'Y0', n_z x n_sec initial biovolume array.
        %> @return result   Struct with fields: time, concentrations, col_grid, profile, cfg.
        function result = run(obj, varargin)

            p = inputParser;
            addParameter(p, 'Y0', [], @isnumeric);
            parse(p, varargin{:});

            Y0 = p.Results.Y0;
            if isempty(Y0)
                Y0 = obj.defaultInitialCondition();
            end

            % time settings
            t0     = obj.cfg.t_init;
            tf     = obj.cfg.t_final;
            dt_out = obj.cfg.delta_t;
            dt     = dt_out;   % internal step = output step for now

            % CFL check
            cfl = ColumnTransport.maxCFL(obj.rhs.w_z, obj.profile.Kz, ...
                                          obj.col_grid.dz, dt);
            if cfl > 0.9
                warning('ColumnSimulation:highCFL', ...
                    'CFL = %.3f > 0.9. Consider smaller dt or coarser dz.', cfl);
            end
            fprintf('ColumnSimulation: CFL = %.4f, dt = %.3f day, dz = %.1f m\n', ...
                cfl, dt, obj.col_grid.dz);

            % output times
            t_out = (t0 : dt_out : tf)';
            n_t   = length(t_out);
            n_z   = obj.col_grid.n_z;
            n_sec = obj.cfg.n_sections;

            % storage: aggregates (n_t x n_z x n_sec), fecal pellets same size
            Y_hist   = zeros(n_t, n_z, n_sec);
            Yfp_hist = zeros(n_t, n_z, n_sec);

            % initial conditions: aggregates from defaultIC, fecal starts at zero
            Y   = Y0;
            Yfp = zeros(n_z, n_sec);

            Y_hist(1, :, :)   = Y;
            Yfp_hist(1, :, :) = Yfp;
            t = t0;
            i_out = 2;

            while t < tf - 1e-10
                [Y, Yfp] = obj.rhs.stepY(Y, dt, Yfp);
                t = t + dt;
                if i_out <= n_t && abs(t - t_out(i_out)) < dt * 0.5
                    Y_hist(i_out, :, :)   = Y;
                    Yfp_hist(i_out, :, :) = Yfp;
                    i_out = i_out + 1;
                end
            end

            % pack result
            result.time                 = t_out;
            result.concentrations       = Y_hist;
            result.fecal_concentrations = Yfp_hist;   % separate fecal pellet array
            result.col_grid             = obj.col_grid;
            result.profile              = obj.profile;
            result.cfg                  = obj.cfg;
            result.cfl                  = cfl;
            result.w_z                  = obj.rhs.w_z;
            result.w_fp_z               = obj.rhs.w_fp_z;
        end

        function Y0 = defaultInitialCondition(obj)
            % Pulse in the top layer: power-law spectrum in layer 1, zeros below.
            n_z   = obj.col_grid.n_z;
            n_sec = obj.cfg.n_sections;
            v0    = InitialSpectrumBuilder.initialSpectrum(obj.cfg, obj.size_grid);
            Y0    = zeros(n_z, n_sec);
            Y0(1, :) = v0(:)';
        end
    end
end
