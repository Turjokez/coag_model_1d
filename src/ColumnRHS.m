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
        w_z         % n_z x n_sec, sinking speed [m/day]

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

            % sinking speed field
            obj.w_z = obj.buildWindField();

            % depth-scaling vectors
            obj.brown_scale = profile.brownianScale(cfg);
            obj.shear_scale = profile.shearScale(cfg);
            obj.ds_scale    = profile.dsScale(cfg);
        end

        function Y_new = stepY(obj, Y, dt)
            % STEPY  One full time step: transport then depth-scaled process rates.
            % Y: n_z x n_sec,  dt: day

            % 1. transport (advection + diffusion)
            Y_new = ColumnTransport.step(Y, obj.w_z, obj.profile.Kz, obj.col_grid.dz, dt);

            % 2. process rates with depth-specific kernel scaling
            % use substeps to keep explicit Euler stable
            n_sub = max(1, round(obj.cfg_orig.proc_substeps));
            dt_sub = dt / n_sub;
            n_z = obj.col_grid.n_z;

            for k = 1:n_z
                sb = obj.brown_scale(k);
                ss = obj.shear_scale(k);
                sd = obj.ds_scale(k);

                % scaled beta matrices for this depth layer (constant across substeps)
                b25_k = sb .* obj.b25_brown + ss .* obj.b25_shear + sd .* obj.b25_ds;
                b1_k  = sb .* obj.b1_brown  + ss .* obj.b1_shear  + sd .* obj.b1_ds;

                v_k = Y_new(k, :)';
                for s = 1:n_sub
                    dvdt = obj.coag_rhs.evaluateScaled(0, v_k, b25_k, b1_k);
                    v_k  = max(v_k + dt_sub * dvdt, 0);
                end
                Y_new(k,:) = v_k';
            end

            % 3. operator_split disagg (if enabled)
            if obj.useOperatorSplitDisagg()
                Y_new = obj.applyDisaggSplit(Y_new);
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
                Y(k,:) = DisaggregationOperatorSplit.apply(v_k, obj.size_grid, obj.cfg_orig)';
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
