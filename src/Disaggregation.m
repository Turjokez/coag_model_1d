classdef Disaggregation
    %DISAGGREGATION Simple turbulence-driven breakup (legacy-style)
    %
    % This implements the same net disaggregation term that is used in the
    % RHS legacy loop. It does NOT enforce mass conservation by itself.
    % For now, we treat this as a physical loss term in diagnostics.
    % The intent here is to keep behavior identical while enabling tests
    % and clean encapsulation.
    %
    % term(isec) = -c3 * c4^isec * (v(isec) - c4*v(isec+1)), isec=2..n-1
    %
    % Inputs:
    %   v   - state vector (n x 1)
    %   cfg - SimulationConfig (uses cfg.c3, cfg.c4)

    methods (Static)
        function term = netTerm(v, cfg)
            n = length(v);
            term = zeros(n,1);

            if n < 3
                return;
            end

            c3 = 0.02;
            c4 = 1.45;
            if nargin >= 2 && ~isempty(cfg)
                if isprop(cfg,'c3') && ~isempty(cfg.c3), c3 = cfg.c3; end
                if isprop(cfg,'c4') && ~isempty(cfg.c4), c4 = cfg.c4; end
            end

            % Clip negatives to zero for physics (no EPS injection)
            v_pos = max(v, 0);

            for isec = 2:(n-1)
                term(isec) = term(isec) - c3 * c4^isec * (v_pos(isec) - c4 * v_pos(isec+1));
            end
        end
    end
end
