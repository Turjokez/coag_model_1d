classdef DisaggregationOperatorSplit
    %DISAGGREGATIONOPERATORSPLIT Operator-split disaggregation (Alldredge-style)
    %
    % Applies a post-step redistribution based on a maximum size (D_max):
    %   - Sum all biovolume in bins >= D_max
    %   - Send disagg_frac_next to the next smaller bin
    %   - Distribute the remainder uniformly across all smaller bins
    %   - Zero out bins >= D_max
    %
    % This matches the logic in legacy ResetInitialCond.m.

    methods (Static)
        function v_out = apply(v_in, grid, cfg, t)
            if nargin < 4
                t = 0;
            end

            v = v_in(:);
            n = length(v);
            if n == 0
                v_out = v_in;
                return;
            end

            % Guard against negative noise
            v = max(v, 0);

            [d_max_idx, ~] = DisaggregationOperatorSplit.maxIndex(grid, cfg, t);

            % Nothing to do if d_max exceeds the largest bin
            if d_max_idx > n
                v_out = v;
                return;
            end

            biovol_disagg = sum(v(d_max_idx:end));

            % If D_max is too small to redistribute in the usual way
            if d_max_idx <= 2
                v(1) = v(1) + biovol_disagg;
                v(d_max_idx:end) = 0.0;
                v_out = v;
                return;
            end

            frac_next = 2/3;
            if isprop(cfg,'disagg_frac_next') && ~isempty(cfg.disagg_frac_next)
                frac_next = cfg.disagg_frac_next;
            end
            frac_next = max(0, min(1, frac_next));

            v(d_max_idx-1) = v(d_max_idx-1) + frac_next * biovol_disagg;

            b_small = (1.0 - frac_next) * biovol_disagg / (d_max_idx-2);
            v(1:d_max_idx-2) = v(1:d_max_idx-2) + b_small;

            v(d_max_idx:end) = 0.0;
            v_out = v;
        end

        function [d_max_idx, d_max_cm] = maxIndex(grid, cfg, t)
            d_max_cm = DisaggregationOperatorSplit.maxDiameterCm(cfg, t);

            d_low = 2.0 * (grid.amfrac * grid.v_lower.^grid.bmfrac);
            d_low_max = max(d_low);
            if d_max_cm >= d_low_max
                d_max_idx = length(d_low) + 1; % no redistribution if Dmax exceeds grid
                return;
            end
            [~, ix] = min(abs(d_low - d_max_cm));
            d_max_idx = ix;
        end

        function d_max_cm = maxDiameterCm(cfg, t)
            d_max_cm = [];
            if isprop(cfg,'disagg_dmax_cm') && ~isempty(cfg.disagg_dmax_cm)
                d_max_cm = cfg.disagg_dmax_cm;
            end

            if isempty(d_max_cm)
                epsilon = DisaggregationOperatorSplit.epsilonAt(cfg, t);
                if isempty(epsilon) || ~isfinite(epsilon)
                    error('Operator-split disagg requires disagg_epsilon or disagg_dmax_cm.');
                end

                c_disagg = 3.0;
                gamma_disagg = 0.15;
                if isprop(cfg,'disagg_C') && ~isempty(cfg.disagg_C)
                    c_disagg = cfg.disagg_C;
                end
                if isprop(cfg,'disagg_gamma') && ~isempty(cfg.disagg_gamma)
                    gamma_disagg = cfg.disagg_gamma;
                end

                % D_max = C * epsilon^(-gamma) in mm, then convert to cm
                d_max_cm = 0.1 * c_disagg * epsilon.^(-gamma_disagg);
            end
        end

        function epsilon = epsilonAt(cfg, t) %#ok<INUSD>
            epsilon = [];
            if isprop(cfg,'disagg_epsilon') && ~isempty(cfg.disagg_epsilon)
                eps_val = cfg.disagg_epsilon;
                if isa(eps_val, 'function_handle')
                    epsilon = eps_val(t);
                else
                    epsilon = eps_val;
                end
            end
        end
    end
end
