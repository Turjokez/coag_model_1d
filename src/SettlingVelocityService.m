classdef SettlingVelocityService < handle
    %SETTLINGVELOCITYSERVICE Service for calculating particle settling velocities
    %
    % Internal model unit: velocity v [cm/s]
    % Paper laws are converted to cm/s inside velocityForSections.

    methods (Static)

        function v = velocity(r, rcons, setcon)
            %VELOCITY Legacy fractal/Stokes-like settling velocity
            % v     = settling velocities [cm/s]
            % r     = particle radii [cm] (fractal radius)
            % rcons = conserved-volume radii [cm]
            % setcon = (2/9)*(delta_rho/rho)*g/kvisc

            v = setcon * rcons .* rcons .* rcons ./ r;
        end

        function v = velocityWithConfig(r, rcons, config, grid)
            %VELOCITYWITHCONFIG Legacy convenience method
            %#ok<INUSD>
            v = SettlingVelocityService.velocity(r, rcons, grid.setcon);
        end

        function v = velocityForSections(grid, config)
            %VELOCITYFORSECTIONS Settling velocity for each size section [cm/s]
            %
            % Supported config.sinking_law:
            %   'current'     : legacy fractal/Stokes-like law (default)
            %   'kriest_8'    : w = 66  * d^0.62  (m/day), d in cm
            %   'kriest_9'    : w = 132 * d^0.62  (m/day), d in cm
            %   'siegel_2025' : w = 20.2* D^0.67  (m/day), D in mm
            %   'kriest_8_capped' : w = min(66*d^0.62, w_max) (m/day), d in cm
            %   'kriest_8_flat'   : Kriest-8 up to D_flat, then constant (m/day)
            %
            % Supported config.sinking_size:
            %   'volume' (default) or 'image'
            %
            % Optional:
            %   config.sinking_scale (dimensionless multiplier; default 1)

            if nargin < 2 || isempty(config)
                config = SimulationConfig();
            end

            if ~isprop(config,'sinking_law') || isempty(config.sinking_law)
                config.sinking_law = 'current';
            end
            if ~isprop(config,'sinking_size') || isempty(config.sinking_size)
                config.sinking_size = 'volume';
            end

            % default scale if missing
            scale = 1.0;
            if isprop(config,'sinking_scale') && ~isempty(config.sinking_scale)
                scale = config.sinking_scale;
            end

            law = lower(string(config.sinking_law));

            switch law
                case "current"
                    r_i = grid.getFractalRadii();     % cm
                    r_v = grid.getConservedRadii();   % cm
                    v   = SettlingVelocityService.velocity(r_i, r_v, grid.setcon); % cm/s

                case "kriest_8"
                    d_cm    = SettlingVelocityService.getDiameterCm(grid, config); % cm
                    w_m_day = 66 .* (d_cm .^ 0.62);   % m/day
                    v       = SettlingVelocityService.mday_to_cms(w_m_day, config);

                case "kriest_9"
                    d_cm    = SettlingVelocityService.getDiameterCm(grid, config); % cm
                    w_m_day = 132 .* (d_cm .^ 0.62);  % m/day
                    v       = SettlingVelocityService.mday_to_cms(w_m_day, config);

                case "kriest_8_capped"
                    d_cm    = SettlingVelocityService.getDiameterCm(grid, config); % cm
                    w_m_day = 66 .* (d_cm .^ 0.62);   % m/day
                    w_max   = 70; % default cap [m/day]
                    if isprop(config,'sinking_w_max_mday') && ~isempty(config.sinking_w_max_mday)
                        w_max = config.sinking_w_max_mday;
                    end
                    w_m_day = min(w_m_day, w_max);
                    v       = SettlingVelocityService.mday_to_cms(w_m_day, config);

                case "kriest_8_flat"
                    d_cm    = SettlingVelocityService.getDiameterCm(grid, config); % cm
                    w_m_day = 66 .* (d_cm .^ 0.62);   % m/day
                    d_flat_cm = 0.1; % default 1 mm
                    if isprop(config,'sinking_d_flat_cm') && ~isempty(config.sinking_d_flat_cm)
                        d_flat_cm = config.sinking_d_flat_cm;
                    end
                    w_flat = 66 .* (d_flat_cm .^ 0.62);
                    w_m_day(d_cm >= d_flat_cm) = w_flat;
                    v       = SettlingVelocityService.mday_to_cms(w_m_day, config);

                case "siegel_2025"
                    d_cm    = SettlingVelocityService.getDiameterCm(grid, config); % cm
                    D_mm    = d_cm * 10;               % mm
                    w_m_day = 20.2 .* (D_mm .^ 0.67);  % m/day
                    v       = SettlingVelocityService.mday_to_cms(w_m_day, config);

                otherwise
                    error("Unknown sinking_law: %s", config.sinking_law);
            end

            % Apply optional multiplier
            v = v * scale;

            % Safety
            v(~isfinite(v)) = 0;
            v(v < 0) = 0;
        end

        % ----------------- helpers -----------------

        function d_cm = getDiameterCm(grid, config)
            %GETDIAMETERCm Return diameter per section [cm] for empirical laws

            use_image = isprop(config,'sinking_size') && strcmpi(config.sinking_size,'image');

            if use_image
                % Needs config.r_to_rg; if missing, fall back to volume diameter
                if isprop(config,'r_to_rg') && ~isempty(config.r_to_rg)
                    d_cm = grid.getImageDiameters(config);
                else
                    d_cm = grid.getVolumeDiameters();
                end
            else
                d_cm = grid.getVolumeDiameters();
            end
        end

        function v_cms = mday_to_cms(w_m_day, config)
            %MDAY_TO_CMS Convert m/day -> cm/s
            v_cms = (w_m_day * 100) / config.day_to_sec;
        end

    end
end
