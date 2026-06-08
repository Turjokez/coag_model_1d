classdef ColumnRHS < handle
    % COLUMNRHS  Right-hand side for the 1-D depth column model.
    %
    % One time step:
    %   1. Transport (upwind advection + flux-form diffusion).
    %   2. Process rates at each depth layer, with depth-scaled kernels.
    %
    % Sinking is handled by transport — coag_rhs has enable_sinking=false.
    %
    % Phase 1: fixed beta matrices, same rates at all depths.
    % Phase 2: brown_scale(k), shear_scale(k), ds_scale(k) applied per layer.
    %
    % Usage:
    %   rhs = ColumnRHS(cfg, size_grid, col_grid, profile);
    %   Y_new = rhs.stepY(Y, dt);

    properties
        cfg         % SimulationConfig (local copy, enable_sinking = false)
        cfg_orig    % original SimulationConfig
        size_grid   % DerivedGrid
        col_grid    % ColumnGrid
        profile     % DepthProfile
        coag_rhs    % CoagulationRHS — process rates without sinking
        w_z         % n_z x n_sec, aggregate sinking speed [m/day]
        w_fp_z      % n_z x n_sec, fecal pellet sinking speed [m/day]

        % depth-scaling vectors (n_z x 1 each)
        brown_scale
        shear_scale
        ds_scale

        % component beta matrices for per-depth scaling
        b25_brown   % n_sec x n_sec
        b1_brown
        b25_shear
        b1_shear
        b25_ds
        b1_ds

        % depth-varying D_max: D_max(k) = A * eps(k)^(-1/4)
        % More turbulence (high eps) -> smaller D_max -> more fragmentation.
        % A calibrated so D_max ~ 1 mm at surface (eps ~ 1e-9 m^2/s^3).
        Dmax_A = 9.39e-6

        zoo        % ZooplanktonGrazing object (empty if disabled)
        cross_coag % FecalCrossCoag object (empty if disabled)
    end

    methods
        function obj = ColumnRHS(cfg, size_grid, col_grid, profile)
            obj.cfg_orig  = cfg;
            obj.size_grid = size_grid;
            obj.col_grid  = col_grid;
            obj.profile   = profile;

            % local config: disable sinking so transport handles it
            obj.cfg                = cfg.copy();
            obj.cfg.enable_sinking = false;
            obj.cfg.box_depth      = [];

            % if operator-split disagg is selected, disable legacy disagg
            % inside CoagulationRHS — same fix as CoagulationSimulation does.
            % Without this, both legacy (c3=0.02) and operator-split run together.
            if isprop(cfg,'enable_disagg') && cfg.enable_disagg && ...
               isprop(cfg,'disagg_mode')   && strcmpi(cfg.disagg_mode,'operator_split')
                obj.cfg.enable_disagg = false;
            end

            % build beta component matrices (used for per-depth scaling)
            % Each component is stored already scaled by its base physical factor
            % (same factors combineAndScale uses) including alpha, but WITHOUT the
            % depth correction.  At depth k the combined matrix is:
            %   b25_k = brown_scale(k)*b25_brown + shear_scale(k)*b25_shear + ds_scale(k)*b25_ds
            % When all scales = 1 this recovers exactly the flat betas from combineAndScale.
            assembler = BetaAssembler(obj.cfg, size_grid);
            if obj.cfg.enable_coag
                b_brown = assembler.computeFor('KernelBrown');
                b_shear = assembler.computeFor('KernelCurSh');
                b_ds    = assembler.computeFor(ColumnRHS.dsKernelName(cfg));
                betas   = assembler.combineAndScale(b_brown, b_shear, b_ds);

                % get the same base scaling factors used inside combineAndScale
                alpha_val = 1.0;
                if isprop(obj.cfg,'alpha') && ~isempty(obj.cfg.alpha)
                    alpha_val = obj.cfg.alpha;
                end
                s_br = 1.0; if isprop(obj.cfg,'scale_brown') && ~isempty(obj.cfg.scale_brown), s_br = obj.cfg.scale_brown; end
                s_sh = 1.0; if isprop(obj.cfg,'scale_shear') && ~isempty(obj.cfg.scale_shear), s_sh = obj.cfg.scale_shear; end
                s_ds = 1.0; if isprop(obj.cfg,'scale_ds')    && ~isempty(obj.cfg.scale_ds),    s_ds = obj.cfg.scale_ds;    end

                f_brown = alpha_val * size_grid.conBr  * obj.cfg.day_to_sec * s_br;
                f_shear = alpha_val * obj.cfg.gamma    * obj.cfg.day_to_sec * s_sh;
                f_ds    = alpha_val * size_grid.setcon * obj.cfg.day_to_sec * s_ds;

                obj.b25_brown = f_brown .* b_brown.b25;
                obj.b1_brown  = f_brown .* b_brown.b1;
                obj.b25_shear = f_shear .* b_shear.b25;
                obj.b1_shear  = f_shear .* b_shear.b1;
                obj.b25_ds    = f_ds    .* b_ds.b25;
                obj.b1_ds     = f_ds    .* b_ds.b1;
            else
                % coag off — build dummy zero components (not used in stepY)
                betas         = assembler.computeFor('KernelBrown');
                obj.b25_brown = zeros(size(betas.b25));
                obj.b1_brown  = zeros(size(betas.b1));
                obj.b25_shear = zeros(size(betas.b25));
                obj.b1_shear  = zeros(size(betas.b1));
                obj.b25_ds    = zeros(size(betas.b25));
                obj.b1_ds     = zeros(size(betas.b1));
            end

            % linear matrix (growth only, no sinking)
            lin    = LinearProcessBuilder.linearMatrix(obj.cfg, size_grid);
            [Dm, Dp] = LinearProcessBuilder.disaggregationMatrices(obj.cfg);
            obj.coag_rhs = CoagulationRHS(betas, lin, Dm, Dp, obj.cfg, size_grid);

            % sinking speed fields
            obj.w_z    = obj.buildWindField();
            obj.w_fp_z = obj.buildFpWindField();

            % depth-scaling vectors
            obj.brown_scale = profile.brownianScale(cfg);
            obj.shear_scale = profile.shearScale(cfg);
            obj.ds_scale    = profile.dsScale(cfg);

            % build zoo object if grazing is enabled
            if isprop(cfg, 'enable_zoo') && cfg.enable_zoo
                obj.zoo = ZooplanktonGrazing( ...
                    'Zc', cfg.zoo_Zc, ...
                    'c',  cfg.zoo_c,  ...
                    'Zf', cfg.zoo_Zf, ...
                    's',  cfg.zoo_s,  ...
                    'p',  cfg.zoo_p,  ...
                    'ic', cfg.zoo_ic);
            else
                obj.zoo = [];
            end

            % build cross-coag object (fecal pellets sticking to marine snow)
            % uses reference sinking speeds at surface (depth scaling applied in stepY)
            if isprop(cfg, 'enable_zoo') && cfg.enable_zoo
                w_fp_ref  = obj.w_fp_z(1, :)';   % surface fecal speed [m/day]
                w_agg_ref = obj.w_z(1, :)';        % surface agg speed  [m/day]
                obj.cross_coag = FecalCrossCoag(cfg, size_grid, w_fp_ref, w_agg_ref);
            else
                obj.cross_coag = [];
            end
        end

        function [Y_new, Yfp_new] = stepY(obj, Y, dt, Yfp)
            % STEPY  One full time step: transport then depth-scaled process rates.
            %
            % Y:    n_z x n_sec  — aggregate biovolume
            % Yfp:  n_z x n_sec  — fecal pellet biovolume (optional, zeros if absent)
            % Both arrays are evolved one step and returned.
            %
            % Fecal pellets (Yfp):
            %   - Undergo the same sinking transport as aggregates (same w_z for now).
            %   - Receive fecal production from zooplankton grazing.
            %   - Do NOT undergo coagulation or disaggregation in this step.

            n_z = obj.col_grid.n_z;

            % initialise Yfp to zeros if not provided (backward compatible)
            if nargin < 4 || isempty(Yfp)
                Yfp = zeros(size(Y));
            end

            % 1. transport aggregates (advection + diffusion)
            Y_new = ColumnTransport.step(Y, obj.w_z, obj.profile.Kz, obj.col_grid.dz, dt);

            % 1b. transport fecal pellets (faster Stokes-based sinking speed)
            Yfp_new = ColumnTransport.step(Yfp, obj.w_fp_z, obj.profile.Kz, obj.col_grid.dz, dt);

            % 2. coagulation process rates at each depth layer
            n_sub  = max(1, round(obj.cfg_orig.proc_substeps));
            dt_sub = dt / n_sub;

            for k = 1:n_z
                sb = obj.brown_scale(k);
                ss = obj.shear_scale(k);
                sd = obj.ds_scale(k);

                b25_k = sb .* obj.b25_brown + ss .* obj.b25_shear + sd .* obj.b25_ds;
                b1_k  = sb .* obj.b1_brown  + ss .* obj.b1_shear  + sd .* obj.b1_ds;

                v_k = Y_new(k, :)';
                for s = 1:n_sub
                    dvdt = obj.coag_rhs.evaluateScaled(0, v_k, b25_k, b1_k);
                    v_k  = max(v_k + dt_sub * dvdt, 0);
                end
                Y_new(k,:) = v_k';
            end

            % 3. operator-split grazing at each depth layer (if enabled)
            % Fecal production goes to Yfp at the correct bin, not back into Y.
            if ~isempty(obj.zoo)
                day_to_sec = obj.cfg_orig.day_to_sec;
                n_sec      = obj.cfg_orig.n_sections;
                target_bin = max(1, min(n_sec, round(obj.zoo.ic) + 1));

                for k = 1:n_z
                    v_k   = Y_new(k, :)';
                    w_cms = obj.w_z(k, :)' .* (100 / day_to_sec);

                    if ~isempty(obj.profile) && ~isempty(obj.profile.Zc)
                        [dvdt, fp_flux] = obj.zoo.graze(v_k, w_cms, ...
                                              obj.profile.Zc(k), obj.profile.Zf(k));
                    else
                        [dvdt, fp_flux] = obj.zoo.graze(v_k, w_cms);
                    end

                    % update aggregate array (losses only, no fecal return here)
                    Y_new(k, :) = max(v_k + dt .* dvdt, 0)';

                    % add fecal production to fecal pellet array at this depth
                    Yfp_new(k, target_bin) = max(0, Yfp_new(k, target_bin) + dt * fp_flux);
                end
            end

            % 3b. cross-coagulation: fecal pellets stick to marine snow at each layer
            % DS dominates (fecal sinks ~17x faster). alpha_cross = 0.5 by default.
            if ~isempty(obj.cross_coag)
                for k = 1:n_z
                    sd = obj.ds_scale(k);
                    [Y_new(k,:), Yfp_new(k,:)] = obj.cross_coag.apply( ...
                        Y_new(k,:)', Yfp_new(k,:)', dt, sd);
                end
            end

            % 3c. micro-zoo mining at each depth layer (if enabled)
            % Particles shrink bin-by-bin; fecal goes to Yfp at target bin.
            if isprop(obj.cfg_orig,'enable_mining') && obj.cfg_orig.enable_mining && ~isempty(obj.zoo)
                av_vol     = obj.size_grid.av_vol(:);   % cm^3 per bin
                day_to_sec = obj.cfg_orig.day_to_sec;
                n_sec      = obj.cfg_orig.n_sections;
                target_bin = max(1, min(n_sec, round(obj.zoo.ic) + 1));

                for k = 1:n_z
                    v_k   = Y_new(k, :)';
                    w_cms = obj.w_z(k, :)' .* (100 / day_to_sec);
                    if ~isempty(obj.profile) && isprop(obj.profile, 'Zm') && ~isempty(obj.profile.Zm)
                        Zm = obj.profile.Zm(k);
                    else
                        Zm = obj.cfg_orig.mining_Zm;
                    end

                    min_bin = obj.cfg_orig.mining_min_bin;
                    [dvdt_m, fp_m] = obj.zoo.mine(v_k, w_cms, av_vol, ...
                                         Zm, ...
                                         obj.cfg_orig.mining_dm, ...
                                         obj.cfg_orig.mining_s, ...
                                         min_bin);

                    Y_new(k, :)            = max(v_k + dt .* dvdt_m, 0)';
                    Yfp_new(k, target_bin) = max(0, Yfp_new(k, target_bin) + dt * fp_m);
                end
            end

            % 4. surface production — layer 1 only, aggregate array only
            if isprop(obj.cfg_orig,'enable_surface_pp') && obj.cfg_orig.enable_surface_pp
                ib = max(1, min(obj.cfg_orig.n_sections, obj.cfg_orig.surface_pp_bin));
                use_mu = isprop(obj.cfg_orig,'surface_pp_mu') && obj.cfg_orig.surface_pp_mu > 0;
                if use_mu
                    Y_new(1, ib) = Y_new(1, ib) * (1 + dt * obj.cfg_orig.surface_pp_mu);
                else
                    Y_new(1, ib) = Y_new(1, ib) + dt * obj.cfg_orig.surface_pp_rate;
                end
            end

            % 5. operator-split disagg on aggregates only
            if obj.useOperatorSplitDisagg()
                Y_new = obj.applyDisaggSplit(Y_new);
            end

            % 6. microbial remineralization (first-order, operator-split)
            % Exact decay Y *= exp(-r*dt) — never produces negatives.
            % r can scale with temperature (Q10) and size (d^-gamma).
            if isprop(obj.cfg_orig,'enable_microbe') && obj.cfg_orig.enable_microbe
                d_cm = obj.size_grid.dcomb(:);   % n_sec x 1, bin diameter [cm]
                for k = 1:n_z
                    % base rate
                    r = obj.cfg_orig.microbe_r0;

                    % optional Q10 temperature scaling
                    if obj.cfg_orig.microbe_use_temp
                        T_C = obj.profile.T_K(k) - 273.15;
                        r = r * obj.cfg_orig.microbe_q10 ^ ...
                            ((T_C - obj.cfg_orig.microbe_tref_C) / 10);
                    end

                    % optional size scaling: smaller particles degrade faster
                    if obj.cfg_orig.microbe_gamma_size ~= 0
                        r_vec = r .* (d_cm / obj.cfg_orig.microbe_dref_cm) .^ ...
                                    (-obj.cfg_orig.microbe_gamma_size);
                    else
                        r_vec = r * ones(obj.cfg_orig.n_sections, 1);
                    end

                    % exact exponential decay (stable for any r*dt)
                    Y_new(k,:)   = Y_new(k,:)'   .* exp(-dt .* r_vec);
                    Yfp_new(k,:) = Yfp_new(k,:)' .* exp(-dt .* r_vec .* obj.cfg_orig.microbe_fp_mult);
                end
            end
        end

        function w = buildWindField(obj)
            % w(k, s): sinking speed [m/day] at depth k for size s.
            % Viscosity correction: w(k) = w_ref * nu_ref / nu(k)
            if isprop(obj.cfg_orig,'enable_sinking') && ~logical(obj.cfg_orig.enable_sinking)
                w = zeros(obj.col_grid.n_z, obj.cfg_orig.n_sections);
                return;
            end
            v_cms  = SettlingVelocityService.velocityForSections(obj.size_grid, obj.cfg_orig);
            w_ref  = (v_cms / 100) * obj.cfg.day_to_sec;  % m/day
            nu_ref = obj.cfg_orig.kvisc;
            scale  = nu_ref ./ obj.profile.nu;             % n_z x 1
            w      = scale * w_ref(:)';                    % n_z x n_sec
        end

        function w = buildFpWindField(obj)
            % Fecal pellet sinking: Stokes law with dense-pellet excess density.
            % Same viscosity correction per depth layer as aggregates.
            if isprop(obj.cfg_orig,'enable_sinking') && ~logical(obj.cfg_orig.enable_sinking)
                w = zeros(obj.col_grid.n_z, obj.cfg_orig.n_sections);
                return;
            end
            v_cms  = SettlingVelocityService.velocityFecalPellets(obj.size_grid, obj.cfg_orig);
            w_ref  = (v_cms / 100) * obj.cfg_orig.day_to_sec;  % m/day
            nu_ref = obj.cfg_orig.kvisc;
            scale  = nu_ref ./ obj.profile.nu;                  % n_z x 1
            w      = scale * w_ref(:)';                         % n_z x n_sec
        end

        function ok = useOperatorSplitDisagg(obj)
            ok = isprop(obj.cfg_orig,'enable_disagg') && obj.cfg_orig.enable_disagg ...
              && isprop(obj.cfg_orig,'disagg_mode') ...
              && strcmpi(obj.cfg_orig.disagg_mode, 'operator_split');
        end

        function Y = applyDisaggSplit(obj, Y)
            % Operator_split disagg applied independently at each depth layer.
            n_z = obj.col_grid.n_z;
            for k = 1:n_z
                v_k    = Y(k, :)';
                cfg_k = obj.cfg_orig;

                % Depth-varying D_max from eps(k) when available.
                if ~isempty(obj.profile) && isprop(obj.profile, 'eps') ...
                        && numel(obj.profile.eps) >= k ...
                        && isfinite(obj.profile.eps(k)) && obj.profile.eps(k) > 0
                    eps_cm = obj.profile.eps(k);        % cm^2/s^3
                    eps_m  = eps_cm / 1e4;              % m^2/s^3
                    dmax_m = obj.Dmax_A * eps_m^(-1/4); % high eps -> small D_max
                    dmax_cm = 100 * dmax_m;

                    cfg_k = obj.cfg_orig.copy();
                    cfg_k.disagg_dmax_cm = dmax_cm;
                end

                Y(k,:) = DisaggregationOperatorSplit.apply(v_k, obj.size_grid, cfg_k)';
            end
        end
    end

    methods (Static)
        function name = dsKernelName(cfg)
            % Pick DS kernel from ds_kernel_mode setting.
            if isprop(cfg, 'ds_kernel_mode') && strcmpi(cfg.ds_kernel_mode, 'sinking_law')
                name = 'KernelCurDSSinkingLaw';
            else
                name = 'KernelCurDS';
            end
        end
    end
end
