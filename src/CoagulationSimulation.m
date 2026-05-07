classdef CoagulationSimulation < handle
    %COAGULATIONSIMULATION Main simulation controller for coagulation system
    properties
        config          % SimulationConfig object
        grid            % DerivedGrid object
        assembler       % BetaAssembler object
        operators       % Struct with linear operators
        rhs             % CoagulationRHS object
        solver          % ODESolver object
        result          % Simulation results struct
    end

    methods
        function obj = CoagulationSimulation(varargin)
            % Constructor:
            %   CoagulationSimulation()
            %   CoagulationSimulation(cfg)
            %   CoagulationSimulation('name',value,...)

            if nargin == 0
                obj.config = SimulationConfig();
            elseif nargin == 1 && isa(varargin{1}, 'SimulationConfig')
                obj.config = varargin{1};
            else
                obj.config = SimulationConfig(varargin{:});
            end

            obj.config.validate();
            obj.initializeComponents();
        end

        function initializeComponents(obj)
            % Build all "static" components that do not depend on time integration

            obj.grid = obj.config.derive();
            obj.assembler = BetaAssembler(obj.config, obj.grid);

            % Build linear operators ONCE (depends on cfg flags inside builder)
            obj.operators = struct();
            obj.operators.growth = LinearProcessBuilder.growthMatrix(obj.config, obj.grid);
            obj.operators.sink_loss = LinearProcessBuilder.sinkingMatrix(obj.config, obj.grid);
            [obj.operators.disagg_minus, obj.operators.disagg_plus] = ...
                LinearProcessBuilder.disaggregationMatrices(obj.config);
            obj.operators.linear = LinearProcessBuilder.linearMatrix(obj.config, obj.grid);

            obj.solver = ODESolver();
            obj.result = struct();
        end

        function result = run(obj, varargin)
            %RUN Execute the coagulation simulation

            fprintf('Starting coagulation simulation...\n');

            % Optional arguments
            p = inputParser;
            addParameter(p, 'tspan', [], @isnumeric);
            addParameter(p, 'v0', [], @isnumeric);
            addParameter(p, 'solver_options', [], @isstruct);
            parse(p, varargin{:});

            % Time span
            if isempty(p.Results.tspan)
                tspan = obj.config.t_init:obj.config.delta_t:obj.config.t_final-1;
            else
                tspan = p.Results.tspan;
            end
            tspan = tspan(:);

            % Initial condition
            if isempty(p.Results.v0)
                v0 = InitialSpectrumBuilder.initialSpectrum(obj.config, obj.grid);
            else
                v0 = p.Results.v0;
            end

            % Beta matrices
            fprintf('Computing coagulation kernels...\n');
            b_brown = obj.assembler.computeFor('KernelBrown');
            b_shear = obj.assembler.computeFor('KernelCurSh');
            b_ds    = obj.assembler.computeFor(obj.dsKernelName());

            betas = obj.assembler.combineAndScale(b_brown, b_shear, b_ds);
            obj.operators.betas = betas; % keep for diagnostics

            % RHS (betas + linear + disagg + config)
            disagg_mode = "legacy";
            if isprop(obj.config,'disagg_mode') && ~isempty(obj.config.disagg_mode)
                disagg_mode = lower(string(obj.config.disagg_mode));
            end
            do_operator_split = false;
            if isprop(obj.config,'enable_disagg') && ~isempty(obj.config.enable_disagg)
                if logical(obj.config.enable_disagg) && disagg_mode == "operator_split"
                    do_operator_split = true;
                end
            end

            rhs_cfg = obj.config;
            if do_operator_split
                rhs_cfg = obj.config.copy();
                rhs_cfg.enable_disagg = false; % prevent legacy disagg in RHS
            end

            obj.rhs = CoagulationRHS( ...
                betas, obj.operators.linear, ...
                obj.operators.disagg_minus, obj.operators.disagg_plus, ...
                rhs_cfg, obj.grid);

            obj.rhs.validate();

            % Solve
            fprintf('Solving ODEs...\n');
            if do_operator_split
                [t, Y] = obj.solveOperatorSplit(obj.rhs, tspan, v0, p.Results.solver_options);
            else
                [t, Y] = obj.solver.solve(obj.rhs, tspan, v0, p.Results.solver_options);
            end

            % Store
            obj.result.time = t;
            obj.result.concentrations = Y;
            obj.result.initial_conditions = v0;
            obj.result.betas = betas;
            obj.result.operators = obj.operators;

            % Diagnostics
            fprintf('Computing diagnostics...\n');
            obj.result.diagnostics = obj.computeDiagnostics(t, Y);

            % Output data
            fprintf('Computing output data...\n');
            obj.result.output_data = OutputGenerator.spectraAndFluxes(t, Y, obj.grid, obj.config);

            fprintf('Simulation completed successfully.\n');
            result = obj.result;
        end

        function diagnostics = computeDiagnostics(obj, t, Y) %#ok<INUSD>
            diagnostics = struct();

            [diagnostics.sectional_gains, diagnostics.sectional_losses] = ...
                MassBalanceAnalyzer.sectional(Y, obj.operators);

            [diagnostics.total_gains, diagnostics.total_losses] = ...
                MassBalanceAnalyzer.total(Y, obj.operators);

            diagnostics.rate_terms = obj.computeRateTermsOverTime(Y);
            diagnostics.mass_conservation = obj.checkMassConservation(Y);
        end

        function rate_terms = computeRateTermsOverTime(obj, Y)
            % NOTE: rhs.rateTerms() MUST return 5 outputs here

            n_times = size(Y, 1);
            n_sec   = size(Y, 2);

            rate_terms = struct();
            rate_terms.term1 = zeros(n_times, n_sec);
            rate_terms.term2 = zeros(n_times, n_sec);
            rate_terms.term3 = zeros(n_times, n_sec);
            rate_terms.term4 = zeros(n_times, n_sec);
            rate_terms.term5 = zeros(n_times, n_sec);

            for i = 1:n_times
                [term1, term2, term3, term4, term5] = obj.rhs.rateTerms(Y(i, :)');
                rate_terms.term1(i, :) = term1';
                rate_terms.term2(i, :) = term2';
                rate_terms.term3(i, :) = term3';
                rate_terms.term4(i, :) = term4';
                rate_terms.term5(i, :) = term5';
            end
        end

        function conservation = checkMassConservation(obj, Y)
            conservation = struct();

            conservation.total_mass = sum(Y, 2);

            if size(Y, 1) > 1
                conservation.mass_change_rate = diff(conservation.total_mass);
                conservation.relative_change = conservation.mass_change_rate ./ conservation.total_mass(1:end-1);
            else
                conservation.mass_change_rate = [];
                conservation.relative_change  = [];
            end

            conservation.is_conserved = all(conservation.total_mass >= 0);

            % Settling loss diagnostic: (Y * diag(sink_loss))
            conservation.settling_loss = [];
            if isfield(obj.operators, 'sink_loss')
                try
                    sink_diag = diag(obj.operators.sink_loss);
                    if isrow(sink_diag), sink_diag = sink_diag'; end

                    if length(sink_diag) == size(Y, 2)
                        conservation.settling_loss = Y * sink_diag;
                    end
                catch ME
                    warning('Could not calculate settling loss: %s', ME.message);
                end
            end
        end

        function kernel_name = dsKernelName(obj)
            mode = "legacy";
            if isprop(obj.config, 'ds_kernel_mode') && ~isempty(obj.config.ds_kernel_mode)
                mode = lower(string(obj.config.ds_kernel_mode));
            end

            if mode == "sinking_law"
                kernel_name = 'KernelCurDSSinkingLaw';
            else
                warning(['ds_kernel_mode is legacy. Using old DS kernel (KernelCurDS). ', ...
                    'Selected sinking_law will not change DS values.']);
                kernel_name = 'KernelCurDS';
            end
        end

        function generateOutputs(obj, plot_flag)
            if nargin < 2, plot_flag = true; end
            if isempty(obj.result), error('No simulation results available. Run simulation first.'); end

            if plot_flag
                fprintf('Generating plots...\n');

                combined_sectional_gains = obj.result.diagnostics.sectional_gains.coag + ...
                    obj.result.diagnostics.sectional_gains.growth;

                combined_sectional_losses = obj.result.diagnostics.sectional_losses.coag + ...
                    obj.result.diagnostics.sectional_losses.settl + ...
                    obj.result.diagnostics.sectional_losses.growth;

                OutputGenerator.plotAll( ...
                    obj.result.time, obj.result.concentrations, ...
                    obj.result.output_data, obj.result.diagnostics.total_gains, ...
                    obj.result.diagnostics.total_losses, obj.config, ...
                    combined_sectional_gains, combined_sectional_losses, ...
                    obj.result.betas);
            end

            obj.displayDiagnosticsSummary();
        end

        function displayDiagnosticsSummary(obj)
            fprintf('\n=== Simulation Diagnostics Summary ===\n');

            fprintf('Simulation time: %.2f to %.2f days\n', obj.result.time(1), obj.result.time(end));
            fprintf('Number of time points: %d\n', length(obj.result.time));
            fprintf('Number of sections: %d\n', size(obj.result.concentrations, 2));

            MassBalanceAnalyzer.displayBalanceSummary( ...
                obj.result.diagnostics.sectional_gains, ...
                obj.result.diagnostics.sectional_losses, ...
                obj.result.time);

            if obj.result.diagnostics.mass_conservation.is_conserved
                fprintf('Mass conservation: PASSED\n');
            else
                fprintf('Mass conservation: FAILED\n');
            end

            if isfield(obj.result, 'betas')
                obj.result.betas.displaySummary();
            end
        end

        function exportResults(obj, filename)
            if nargin < 2
                filename = sprintf('coagulation_simulation_%s.mat', datestr(now, 'yyyymmdd_HHMMSS'));
            end

            OutputGenerator.exportData(obj.result.output_data, filename);

            results_filename = strrep(filename, '.mat', '_full.mat');
            save(results_filename, 'obj');
            fprintf('Full simulation results saved to: %s\n', results_filename);
        end

        function enableTracer(obj)
            warning('Tracer integration not yet implemented');
        end

        function [t_out, Y_out] = solveOperatorSplit(obj, rhs, tspan, v0, solver_options)
            % Operator-split disaggregation: integrate without disagg,
            % then redistribute above D_max after each outer step.

            if length(tspan) < 1
                error('tspan must have at least one time value');
            end
            if any(diff(tspan) <= 0)
                error('tspan must be strictly increasing for operator_split');
            end

            t_out = tspan(:);
            n_out = length(t_out);
            n_sections = length(v0);
            Y_out = zeros(n_out, n_sections);

            current_t = t_out(1);
            current_v = v0(:);
            Y_out(1,:) = current_v';

            % Determine outer timestep
            outer_dt = obj.config.delta_t;
            if isprop(obj.config,'disagg_outer_dt') && ~isempty(obj.config.disagg_outer_dt)
                outer_dt = obj.config.disagg_outer_dt;
            end
            if outer_dt <= 0
                error('disagg_outer_dt must be > 0');
            end

            tol = 1e-10;
            if n_out > 1
                output_dt = t_out(2) - t_out(1);
                ratio = output_dt / outer_dt;
                if abs(ratio - round(ratio)) > 1e-8
                    error('operator_split requires output dt to be an integer multiple of disagg_outer_dt');
                end
                if outer_dt > output_dt + tol
                    error('disagg_outer_dt must be <= output dt');
                end
            end

            next_out = 2;
            t_end = t_out(end);

            while current_t < t_end - tol
                t_next = min(current_t + outer_dt, t_end);

                [t_step, y_step] = obj.solver.solve(rhs, [current_t t_next], current_v, solver_options);
                current_t = t_step(end);
                current_v = y_step(end,:)';

                % Apply operator-split redistribution
                current_v = DisaggregationOperatorSplit.apply(current_v, obj.grid, obj.config, current_t);

                if next_out <= n_out && abs(current_t - t_out(next_out)) <= tol
                    Y_out(next_out,:) = current_v';
                    next_out = next_out + 1;
                end
            end

            if next_out <= n_out
                Y_out(next_out:end,:) = repmat(current_v', n_out - next_out + 1, 1);
            end
        end
    end
end
