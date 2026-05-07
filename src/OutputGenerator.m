classdef OutputGenerator < handle
    %OUTPUTGENERATOR Generates outputs and visualizations for coagulation simulation

    methods (Static)
        function output_data = spectraAndFluxes(t, Y, grid, config)
            %SPECTRAANDFLUXES Compute particle spectra and fluxes
            % t = time vector
            % Y = concentration matrix
            % grid = DerivedGrid object
            % config = SimulationConfig object
            % Returns: struct with computed data

            n_times    = length(t);
            n_sections = length(grid.v_lower);

            % Initialize output arrays (legacy)
            nspec_v  = zeros(n_times, n_sections);
            masspec_v = nspec_v;
            fluxsect  = nspec_v;
            fluxspec  = nspec_v;

            % NEW (0-D bookkeeping): sinking loss rates in RHS units
            sinkLossSect  = nspec_v;            % state units/day per bin
            sinkLossTotal = zeros(n_times, 1);  % state units/day total

            % Get radii and diameters (legacy)
            r_i    = grid.getFractalRadii(); %#ok<NASGU>
            r_v    = grid.getConservedRadii(); %#ok<NASGU>
            diam_i = grid.getImageDiameters(config);
            diam_v = grid.getVolumeDiameters();

            % Settling velocities (legacy)
            v_cms   = SettlingVelocityService.velocityForSections(grid, config); % cm/s
            if isprop(config,'enable_sinking') && ~config.enable_sinking
                v_cms = zeros(size(v_cms));
            end
            set_vel = (v_cms / 100) * config.day_to_sec; % m/day

            % NEW: depth scale for 0-D sinking bookkeeping
            H = config.dz;
            if isprop(config,'box_depth') && ~isempty(config.box_depth)
                H = config.box_depth;
            end
            H = max(H, eps);
            lambda = set_vel / H;  % 1/day (vector)

            % Create time matrices for vectorized ops (legacy)
            diam_i = diam_i';
            diam_v = diam_v';
            diam_i_mat = diam_i(ones(n_times, 1), :);
            diam_v_mat = diam_v(ones(n_times, 1), :);

            % Compute spectra and fluxes for each time (legacy + new sink loss)
            for jindx = 1:n_times
                yout = Y(jindx, :);

                nspec_v(jindx, :)   = yout ./ (1.5 * grid.v_lower') ./ grid.dwidth';
                masspec_v(jindx, :) = yout ./ grid.dwidth';
                fluxsect(jindx, :)  = yout .* set_vel' * 1e6;
                fluxspec(jindx, :)  = masspec_v(jindx, :) .* set_vel' * 1e6;

                % NEW: RHS-consistent sinking loss (state/day)
                sinkLossSect(jindx, :)  = yout .* lambda';
                sinkLossTotal(jindx, 1) = sum(sinkLossSect(jindx, :));
            end

            % Total quantities (legacy)
            total_flux = sum(fluxsect, 2);
            total_mass = sum(Y, 2);

            % Image-based spectra (legacy)
            diaratio   = (config.fr_dim/3) * diam_v_mat ./ diam_i_mat;
            nspec_i    = nspec_v .* diaratio;
            masspec_i  = masspec_v .* diaratio;
            fluxspec_i = fluxspec .* diaratio;

            % Assemble output data (KEEP CORE OUTPUTS OUTSIDE try/catch)
            output_data = struct();
            output_data.nspec_v     = nspec_v;
            output_data.masspec_v   = masspec_v;
            output_data.fluxsect    = fluxsect;
            output_data.fluxspec    = fluxspec;
            output_data.total_flux  = total_flux;
            output_data.total_mass  = total_mass;

            output_data.diam_i      = diam_i;
            output_data.diam_v      = diam_v;
            output_data.diam_i_mat  = diam_i_mat;
            output_data.diam_v_mat  = diam_v_mat;

            output_data.nspec_i     = nspec_i;
            output_data.masspec_i   = masspec_i;
            output_data.fluxspec_i  = fluxspec_i;

            output_data.set_vel     = set_vel;
            output_data.v_lower     = grid.v_lower;
            output_data.dwidth      = grid.dwidth;
            output_data.diaratio    = diaratio;
            output_data.fr_dim      = config.fr_dim;

            % NEW: RHS-consistent sinking diagnostics (0-D)
            output_data.sinkLossSect  = sinkLossSect;
            output_data.sinkLossTotal = sinkLossTotal;
            output_data.lambda        = lambda;
            output_data.H_used        = H;

            % Optional diagnostics snapshot (do not put core outputs here)
            try
                diag = struct();
                diag.t = t;
                diag.Y = Y;
                diag.diam_i = diam_i;
                diag.diam_v = diam_v;
                diag.set_vel = set_vel;
                diag.v_lower = grid.v_lower;
                diag.dwidth  = grid.dwidth;
                diag.nspec_v = nspec_v;
                diag.masspec_v = masspec_v;
                diag.fluxspec = fluxspec;
                diag.diaratio = diaratio;
                diag.fluxspec_i = fluxspec_i;
                save('plot_diag_oop.mat','diag');
            catch
            end
        end

        function plotAll(t, Y, output_data, gains, losses, config, sectional_gains, sectional_losses, betas)
            %PLOTALL Generate all standard plots
            OutputGenerator.plotSpectraAndFluxes(t, Y, output_data);
            OutputGenerator.plotMassBalance(t, gains, losses);
            OutputGenerator.plotCoagVsSettRatio(t, losses);
            OutputGenerator.plot3DSurfaces(t, Y, sectional_gains, sectional_losses, output_data);
        end

        function plotCoagVsSettRatio(t, losses)
            figure(3);
            if isfield(losses, 'coag') && isfield(losses, 'sett')
                ratio = losses.coag ./ max(losses.sett, eps);
                plot(t, ratio);
                set(gca, 'FontName', 'Helvetica', 'FontSize', 14)
                xlabel('Time [d]', 'FontName', 'Helvetica', 'FontSize', 14)
                ylabel('(Coag Losses)/(Settling Losses)', 'FontName', 'Helvetica', 'FontSize', 14)
            end
        end

        function plotSpectraAndFluxes(t, Y, output_data)
            figure(1);

            subplot(2, 2, 1);
            ns_init  = output_data.nspec_i(1, :);
            ns_final = output_data.nspec_i(end, :);
            loglog(output_data.diam_i, ns_init, 'b', output_data.diam_i, ns_final, 'r');
            xlabel('Particle diameter [cm]');
            ylabel('Number spectrum [# cm^{-4}]');
            axis tight;

            subplot(2, 2, 3);
            plot(output_data.diam_i_mat(:,2:end)', output_data.fluxspec_i(:,2:end)');
            xlabel('Particle image diameter [cm]');
            ylabel('Volume flux spectra [cm^2 m^{-2} d^{-1}]');
            axis tight;

            subplot(2, 4, 3);
            semilogy(t, Y, t, output_data.total_mass, '*--');
            xlabel('Time [d]');
            ylabel('Sectional concentration [vol/vol/sect]');

            subplot(2, 4, 7);
            plot(t, output_data.fluxsect, t, output_data.total_flux, '*--');
            xlabel('Time [d]');
            ylabel('Sectional Flux [cm^3 m^{-2} d^{-1} sect^{-1}]');

            subplot(2, 4, 8);
            plot(t, output_data.total_flux ./ output_data.total_mass / 1e6);
            xlabel('Time [d]');
            ylabel('Average v [m d^{-1}]');
        end

        function plotMassBalance(t, gains, losses)
            figure(2);
            clf;
            if isfield(gains, 'growth') && isfield(losses, 'sett') && isfield(losses, 'coag')
                denominator = losses.sett + losses.coag;
                ratio = gains.growth ./ max(denominator, eps);
                plot(t, ratio);
                xlabel('Time [d]');
                ylabel('Gains/Losses');
                title('Total System Mass Balance');
            end
        end

        function plot3DSurfaces(t, Y, gains, losses, output_data)
            n_sections = size(Y, 2);
            if length(t) ~= size(Y, 1) || n_sections == 0
                warning('Dimension mismatch in plot3DSurfaces: t=%d, Y=%dx%d', length(t), size(Y, 1), n_sections);
                return;
            end

            t_mat = t(:, ones(1, n_sections));

            if length(output_data.v_lower) == n_sections
                v_mat = output_data.v_lower';
                v_mat = v_mat(ones(length(t), 1), :);
            else
                warning('Dimension mismatch: v_lower length=%d, n_sections=%d', length(output_data.v_lower), n_sections);
                v_mat = ones(length(t), n_sections);
            end

            if size(gains, 1) == length(t) && size(gains, 2) == n_sections && ...
               size(losses, 1) == length(t) && size(losses, 2) == n_sections

                ratio = gains ./ max(losses, eps);
                ratio(losses <= 0) = NaN;
                ratio(~isfinite(ratio)) = NaN;
                ratio(ratio < 0) = NaN;

                figure(4);
                surf(log10(v_mat), t_mat, ratio);
                xlabel('Volume');
                ylabel('Time [d]');
                zlabel('Gains/Losses');

                valid_ratios = ratio(~isnan(ratio) & ratio > 0);
                if ~isempty(valid_ratios)
                    zlim([0, max(valid_ratios)]);
                end
            end
        end

        function exportData(output_data, filename)
            if nargin < 2
                filename = sprintf('coagulation_output_%s.mat', datestr(now, 'yyyymmdd_HHMMSS'));
            end
            save(filename, 'output_data');
            fprintf('Output data saved to: %s\n', filename);
        end
    end
end
