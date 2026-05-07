classdef CoagulationRHS < handle
    % COAGULATIONRHS ODE right-hand side for coagulation equations
    properties
        betas;          % BetaMatrices object, contains coagulation parameters
        linear;         % Linear matrix (growth - sinking), units 1/day
        disaggMinus;    % Disaggregation loss matrix
        disaggPlus;     % Disaggregation gain matrix
        config;         % SimulationConfig object
        grid;           % DerivedGrid (optional)
    end

    methods
        function obj = CoagulationRHS(betas, linear, disaggMinus, disaggPlus, config, varargin)
            obj.betas       = betas;
            obj.linear      = linear;
            obj.disaggMinus = disaggMinus;
            obj.disaggPlus  = disaggPlus;
            obj.config      = config;
            obj.grid        = [];
            if nargin >= 6 && ~isempty(varargin)
                obj.grid = varargin{1};
            end
        end

        function dvdt = evaluate(obj, t, v) %#ok<INUSD>
            n_sections = length(v);

            % Boundary case: allow PP to inject even if all zeros
            if all(v == 0)
                dvdt = zeros(size(v));
                if isprop(obj.config,'enable_pp') && obj.config.enable_pp
                    source = 0;
                    if isprop(obj.config,'pp_source') && ~isempty(obj.config.pp_source) && obj.config.pp_source ~= 0
                        source = obj.config.pp_source;
                    elseif isprop(obj.config,'pp_rate') && ~isempty(obj.config.pp_rate) && obj.config.pp_rate ~= 0
                        source = obj.config.pp_rate;
                    end
                    if source ~= 0
                        ib = 1;
                        if isprop(obj.config,'pp_bin') && ~isempty(obj.config.pp_bin)
                            ib = obj.config.pp_bin;
                        end
                        ib = max(1, min(n_sections, ib));
                        dvdt(ib) = dvdt(ib) + source;
                    end
                end
                return
            end

            % numerical guard for internal math
            v_pos = max(v, eps);

            v_r     = v_pos.';                 % row
            v_shift = [0, v_r(1:n_sections-1)]; % row

            % ---- switches (with safe defaults) ----
            do_coag = true;
            if isprop(obj.config,'enable_coag') && ~isempty(obj.config.enable_coag)
                do_coag = logical(obj.config.enable_coag);
            end

            do_lin = true;
            if isprop(obj.config,'enable_linear') && ~isempty(obj.config.enable_linear)
                do_lin = logical(obj.config.enable_linear);
            end
            % If sinking is enabled, keep linear term on even if enable_linear is false
            if isprop(obj.config,'enable_sinking') && ~isempty(obj.config.enable_sinking) && logical(obj.config.enable_sinking)
                do_lin = true;
            end

            do_dis = false;
            if isprop(obj.config,'enable_disagg') && ~isempty(obj.config.enable_disagg)
                do_dis = logical(obj.config.enable_disagg);
            end

            % ---- initialize contributions ----
            term_coag = zeros(n_sections,1);
            term_lin  = zeros(n_sections,1);
            term_dis  = zeros(n_sections,1);

            % ---- coagulation ----
            if do_coag
                % term1 = v_i * sum_j (beta25_ij * v_j)
                term1 = v_r * obj.betas.b25;  % row
                term1 = v_r .* term1;         % row

                % term2 = v_{i-1} * sum_j (b1_{i-1,j} * v_j) (implemented as legacy)
                term2 = v_r * obj.betas.b1;   % row
                term2 = term2 .* v_shift;     % row

                term_coag = (term1 + term2).'; % column
            end

            % ---- linear (growth - sinking) ----
            if do_lin
                term_lin = obj.linear * v_pos; % column
            end

            % ---- disaggregation ----
            if do_dis
                term_dis = Disaggregation.netTerm(v, obj.config);
            end

            % ---- combine ----
            dvdt = term_coag + term_lin + term_dis;

            % ---- PP source (constant injection into one bin) ----
            if isprop(obj.config,'enable_pp') && obj.config.enable_pp
                source = 0;
                if isprop(obj.config,'pp_source') && ~isempty(obj.config.pp_source) && obj.config.pp_source ~= 0
                    source = obj.config.pp_source;
                elseif isprop(obj.config,'pp_rate') && ~isempty(obj.config.pp_rate) && obj.config.pp_rate ~= 0
                    source = obj.config.pp_rate;
                end

                if source ~= 0
                    ib = 1;
                    if isprop(obj.config,'pp_bin') && ~isempty(obj.config.pp_bin)
                        ib = obj.config.pp_bin;
                    end
                    ib = max(1, min(n_sections, ib));
                    dvdt(ib) = dvdt(ib) + source;
                end
            end
        end

        function dvdt = evaluateScaled(obj, t, v, b25_scaled, b1_scaled) %#ok<INUSD>
            % Like evaluate(), but uses caller-supplied beta matrices.
            % Used by ColumnRHS to apply per-depth kernel scaling.
            % b25_scaled, b1_scaled: depth-scaled combined matrices (n_sec x n_sec).
            n_sections = length(v);

            if all(v == 0)
                dvdt = zeros(size(v));
                return
            end

            v_pos   = max(v, eps);
            v_r     = v_pos.';
            v_shift = [0, v_r(1:n_sections-1)];

            do_coag = true;
            if isprop(obj.config,'enable_coag') && ~isempty(obj.config.enable_coag)
                do_coag = logical(obj.config.enable_coag);
            end
            do_lin = true;
            if isprop(obj.config,'enable_linear') && ~isempty(obj.config.enable_linear)
                do_lin = logical(obj.config.enable_linear);
            end
            do_dis = false;
            if isprop(obj.config,'enable_disagg') && ~isempty(obj.config.enable_disagg)
                do_dis = logical(obj.config.enable_disagg);
            end

            term_coag = zeros(n_sections, 1);
            term_lin  = zeros(n_sections, 1);
            term_dis  = zeros(n_sections, 1);

            if do_coag
                term1     = v_r .* (v_r * b25_scaled);
                term2     = (v_r * b1_scaled) .* v_shift;
                term_coag = (term1 + term2).';
            end
            if do_lin
                term_lin = obj.linear * v_pos;
            end
            if do_dis
                term_dis = Disaggregation.netTerm(v, obj.config);
            end

            dvdt = term_coag + term_lin + term_dis;

            % PP source (same behavior as evaluate)
            if isprop(obj.config,'enable_pp') && obj.config.enable_pp
                source = 0;
                if isprop(obj.config,'pp_source') && ~isempty(obj.config.pp_source) && obj.config.pp_source ~= 0
                    source = obj.config.pp_source;
                elseif isprop(obj.config,'pp_rate') && ~isempty(obj.config.pp_rate) && obj.config.pp_rate ~= 0
                    source = obj.config.pp_rate;
                end
                if source ~= 0
                    ib = 1;
                    if isprop(obj.config,'pp_bin') && ~isempty(obj.config.pp_bin)
                        ib = obj.config.pp_bin;
                    end
                    ib = max(1, min(n_sections, ib));
                    dvdt(ib) = dvdt(ib) + source;
                end
            end
        end

        function J = jacobian(obj, t, v) %#ok<INUSD>
            % For now, keep your existing analytic Jacobian but apply switches.
            % This avoids mismatch if coag/linear/disagg are disabled.

            n_sections = length(v);
            v_r = v.'; % row

            % defaults
            do_coag = true;
            if isprop(obj.config,'enable_coag') && ~isempty(obj.config.enable_coag)
                do_coag = logical(obj.config.enable_coag);
            end

            do_lin = true;
            if isprop(obj.config,'enable_linear') && ~isempty(obj.config.enable_linear)
                do_lin = logical(obj.config.enable_linear);
            end
            if isprop(obj.config,'enable_sinking') && ~isempty(obj.config.enable_sinking) && logical(obj.config.enable_sinking)
                do_lin = true;
            end

            do_dis = false;
            if isprop(obj.config,'enable_disagg') && ~isempty(obj.config.enable_disagg)
                do_dis = logical(obj.config.enable_disagg);
            end

            % start with zeros
            J = zeros(n_sections);

            if do_coag
                v_mat   = v_r(ones(1, n_sections), :);
                v_shift = [zeros(n_sections, 1), v_mat(:, 1:end-1)];

                term1 = v_r * obj.betas.b25;
                term1 = diag(term1) + diag(v_r) .* obj.betas.b25;

                term2a = v_r * obj.betas.b1;
                term2a = diag(term2a(2:end), -1);

                term2b = diag(obj.betas.b1, 1);
                term2b = term2b' .* v_r(1:end-1);
                term2b = diag(term2b, -1);

                term2c = diag(v_r(2:end), -1) .* obj.betas.b25';
                term2  = term2a + term2b + term2c;

                term3a = obj.betas.b1  .* v_shift;
                term3b = obj.betas.b25 .* v_mat;
                term3  = (term3a + term3b)';
                term3  = triu(term3, 2) + tril(term3, -1);

                J = J + term1 + term2 + term3;
            end

            if do_lin
                J = J + obj.linear;
            end

            % If you keep matrix disagg in the future, apply here.
            % Your legacy loop is not included in Jacobian (it wasn't before in a clean way).
            if do_dis
                % Optional: if you want matrix-form disagg, you can add:
                % J = J - obj.disaggMinus + obj.disaggPlus;
                % But only do this if your evaluate() uses these matrices too.
            end

            % PP source has zero derivative w.r.t state (constant), so no Jacobian change.
        end

        function [term1, term2, term3, term4, term5] = rateTerms(obj, v)
            % RATETERMS (legacy-compatible):
            % term1 = coag part A (b25)
            % term2 = coag part B (b1 * shift)
            % term3 = linear term (growth - sinking)
            % term4 = disagg net term (legacy loop; if disabled -> 0)
            % term5 = PP source term (if disabled -> 0)
            %
            % All returned as column vectors (n_sections x 1)
        
            n = length(v);
        
            term1 = zeros(n,1);
            term2 = zeros(n,1);
            term3 = zeros(n,1);
            term4 = zeros(n,1);
            term5 = zeros(n,1);
        
            if all(v == 0)
                return
            end
        
            v_pos = max(v, eps);
            v_r     = v_pos.';               % row
            v_shift = [0, v_r(1:n-1)];       % row
        
            % switches
            do_coag = true;
            if isprop(obj.config,'enable_coag') && ~isempty(obj.config.enable_coag)
                do_coag = logical(obj.config.enable_coag);
            end
        
            do_lin = true;
            if isprop(obj.config,'enable_linear') && ~isempty(obj.config.enable_linear)
                do_lin = logical(obj.config.enable_linear);
            end
            if isprop(obj.config,'enable_sinking') && ~isempty(obj.config.enable_sinking) && logical(obj.config.enable_sinking)
                do_lin = true;
            end
        
            do_dis = false;
            if isprop(obj.config,'enable_disagg') && ~isempty(obj.config.enable_disagg)
                do_dis = logical(obj.config.enable_disagg);
            end
        
            % coagulation split into 2 legacy pieces
            if do_coag
                t1 = v_r .* (v_r * obj.betas.b25);     % row
                t2 = (v_r * obj.betas.b1) .* v_shift;  % row
                term1 = t1.';                          % col
                term2 = t2.';                          % col
            end
        
            % linear (growth - sinking)
            if do_lin
                term3 = obj.linear * v_pos;
            end
        
            % disaggregation (net effect)
            if do_dis
                term4 = Disaggregation.netTerm(v, obj.config);
            end
        
            % PP source (constant injection)
            if isprop(obj.config,'enable_pp') && obj.config.enable_pp
                source = 0;
                if isprop(obj.config,'pp_source') && ~isempty(obj.config.pp_source) && obj.config.pp_source ~= 0
                    source = obj.config.pp_source;
                elseif isprop(obj.config,'pp_rate') && ~isempty(obj.config.pp_rate) && obj.config.pp_rate ~= 0
                    source = obj.config.pp_rate;
                end

                if source ~= 0
                    ib = 1;
                    if isprop(obj.config,'pp_bin') && ~isempty(obj.config.pp_bin)
                        ib = obj.config.pp_bin;
                    end
                    ib = max(1, min(n, ib));
                    term5(ib) = term5(ib) + source;
                end
            end
        end

        function validate(obj)
            if isempty(obj.betas) || isempty(obj.linear)
                error('RHS not properly initialized');
            end

            n_sections = obj.config.n_sections;
            if size(obj.linear, 1) ~= n_sections || size(obj.linear,2) ~= n_sections
                error('Linear matrix dimension mismatch');
            end

            if obj.betas.getNumSections() ~= n_sections
                error('Beta matrices dimension mismatch');
            end
        end
    end
end
