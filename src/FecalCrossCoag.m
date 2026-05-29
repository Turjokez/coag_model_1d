classdef FecalCrossCoag
    % FECALCROSSCOAG  Cross-coagulation: fecal pellets stick to marine snow.
    %
    % Rule (Jokulsdottir 2011):
    %   fecal pellet + marine snow -> marine snow (fecal absorbed into aggregate)
    %   fecal pellets do NOT coagulate with each other
    %
    % Physics:
    %   Differential settling (DS) dominates because fecal sinks ~17x faster.
    %   Kernel: curvilinear DS using |w_fp_i - w_agg_j| for the velocity term.
    %   Stickiness: alpha_cross (tunable, default 0.5 -- fecal is compact/low-TEP).
    %
    % Volume budget per collision event (fp bin i + agg bin j):
    %   Y_fp(i)  loses v_i   (fecal particle consumed)
    %   Y(j)     loses v_j   (marine snow particle consumed)
    %   Y(k)     gains v_i+v_j  (new merged marine snow at bin k)
    %   => total volume is conserved
    %
    % Depth scaling:
    %   DS beta scales with ds_scale(k) = nu_ref/nu(k), same as regular DS.

    properties
        n_sec
        alpha_cross    % fecal-marine snow stickiness

        % precomputed matrices at reference (surface) conditions [1/day]
        % dimension: n_sec x n_sec  (row=fp_bin_i, col=agg_bin_j)
        B_loss_fp      % alpha * beta_cross(i,j) / v(j)
        B_loss_agg     % alpha * beta_cross(i,j) / v(i)
        B_gain         % alpha * beta_cross(i,j) * (v(i)+v(j)) / (v(i)*v(j))

        % target marine snow bin when fp_i + agg_j merge (n_sec x n_sec, int)
        % value = n_sec+1 means result overflows grid (lost from top bin)
        target_bin
    end

    methods
        function obj = FecalCrossCoag(cfg, size_grid, w_fp_ref_mday, w_agg_ref_mday)
            % w_fp_ref_mday:  n_sec x 1, fecal sinking speed [m/day] at surface
            % w_agg_ref_mday: n_sec x 1, marine snow sinking speed [m/day] at surface

            n = cfg.n_sections;
            obj.n_sec = n;

            obj.alpha_cross = 0.5;
            if isprop(cfg,'fp_alpha_cross') && ~isempty(cfg.fp_alpha_cross)
                obj.alpha_cross = cfg.fp_alpha_cross;
            end

            v   = size_grid.av_vol(:);     % average bin volume [cm^3], n x 1
            r   = size_grid.dcomb(:) / 2;  % bin radius [cm],           n x 1

            % --- cross-DS beta [cm^3/day] at reference conditions ---
            % curvilinear DS: beta(i,j) = 0.5 * pi * r_min(i,j)^2 * |dw_cmday|
            dw_cmday = abs(w_fp_ref_mday(:) - w_agg_ref_mday(:)') * 100; % cm/day

            r_min = min(r * ones(1,n), ones(n,1) * r');   % n x n
            beta_ref = obj.alpha_cross .* 0.5 * pi .* r_min.^2 .* dw_cmday;

            % --- volume-weighted loss/gain matrices ---
            % divide by v to go from number rate to volume rate
            v_row = v * ones(1,n);   % v(i) repeated across columns
            v_col = ones(n,1) * v';  % v(j) repeated across rows

            obj.B_loss_fp  = beta_ref ./ v_col;              % loss rate of Y_fp(i)
            obj.B_loss_agg = beta_ref ./ v_row;              % loss rate of Y(j)
            obj.B_gain     = beta_ref .* (v_row + v_col) ./ (v_row .* v_col);

            % --- precompute target bins ---
            v_lower = size_grid.v_lower(:);
            v_upper = size_grid.v_upper(:);
            tb = zeros(n, n, 'int32');
            for i = 1:n
                for j = 1:n
                    v_sum = v(i) + v(j);
                    k = find(v_sum >= v_lower & v_sum < v_upper, 1);
                    if isempty(k)
                        tb(i,j) = int32(n+1);  % overflows grid
                    else
                        tb(i,j) = int32(k);
                    end
                end
            end
            obj.target_bin = tb;
        end

        function [Y_new, Yfp_new, vol_xfer] = apply(obj, Y, Yfp, dt, ds_scale)
            % Apply cross-coagulation for one depth layer, one time step.
            %
            % Inputs:
            %   Y, Yfp:   n_sec x 1 volume concentrations
            %   dt:       time step [day]
            %   ds_scale: depth scaling for DS (nu_ref/nu at this layer)
            %
            % Outputs:
            %   Y_new, Yfp_new: updated concentrations
            %   vol_xfer:       total volume moved from Y_fp to Y [same units as Y]

            n = obj.n_sec;
            Y   = Y(:);
            Yfp = Yfp(:);

            % scale precomputed matrices by ds_scale
            Bfp  = ds_scale .* obj.B_loss_fp;
            Bagg = ds_scale .* obj.B_loss_agg;
            Bg   = ds_scale .* obj.B_gain;

            % rate matrix R(i,j): volume transferred per day from pair (i,j)
            outer = Yfp * Y';   % n x n outer product

            % loss from fecal at each bin i: sum over all j
            dYfp = -dt .* sum(Bfp .* outer, 2);   % n x 1

            % loss from marine snow at each bin j: sum over all i
            dY_loss = -dt .* sum(Bagg .* outer, 1)';  % n x 1

            % gain in marine snow at target bins
            R_gain = dt .* Bg .* outer;  % n x n
            tb     = double(obj.target_bin(:));
            valid  = tb <= n;
            dY_gain = accumarray(tb(valid), R_gain(valid), [n, 1]);

            % apply (clamp to zero)
            Yfp_new = max(Yfp + dYfp, 0);
            Y_new   = max(Y + dY_loss + dY_gain, 0);

            % total volume transferred (for budget check)
            vol_xfer = -sum(dYfp);
        end
    end
end
