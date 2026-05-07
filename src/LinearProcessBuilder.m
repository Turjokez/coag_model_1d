classdef LinearProcessBuilder < handle
    methods (Static)

        function G = growthMatrix(cfg, grid) %#ok<INUSD>
            n = cfg.n_sections;

            growth_loss = zeros(n,1);
            growth_gain = zeros(n-1,1);

            if cfg.gro_sec > 0
                growth_loss(cfg.gro_sec:n-1) = -1;
                growth_gain(cfg.gro_sec:end) = 2;
            end

            G = diag(growth_loss) + diag(growth_gain, -1);
            G(1,1) = 1;
            G = cfg.growth * G;
        end

        function S = sinkingMatrix(cfg, grid)
            n = cfg.n_sections;
        
            % hard OFF switch
            if isprop(cfg,'enable_sinking') && ~cfg.enable_sinking
                S = zeros(n);
                return
            end
        
            % v is cm/s
            v_cms   = SettlingVelocityService.velocityForSections(grid, cfg);
        
            % convert to m/day (same as OutputGenerator)
            w_mday  = (v_cms/100) * cfg.day_to_sec;
        
            H = cfg.dz;
            if isprop(cfg,'box_depth') && ~isempty(cfg.box_depth)
                H = cfg.box_depth;
            end
            H = max(H, eps);
        
            settling_rate = w_mday / H;   % 1/day
            S = diag(settling_rate);
        end

        function [Dminus, Dplus] = disaggregationMatrices(cfg)
            n = cfg.n_sections;
            Dminus = zeros(n);
            Dplus  = zeros(n);

            if n > 2
                for k = 2:(n-1)
                    Dminus(k,k)  = cfg.c3 * cfg.c4^k;
                    Dplus(k,k-1) = cfg.c3 * cfg.c4^(k+1);
                end
            end
        end

        function L = linearMatrix(cfg, grid)
            n = cfg.n_sections;
            L = zeros(n);

            % growth part (controlled by enable_linear)
            if ~isprop(cfg,'enable_linear') || cfg.enable_linear
                if isprop(cfg,'growth') && cfg.growth ~= 0
                    L = L + LinearProcessBuilder.growthMatrix(cfg, grid);
                end
            end

            % sinking part (controlled by enable_sinking)
            if ~isprop(cfg,'enable_sinking') || cfg.enable_sinking
                L = L - LinearProcessBuilder.sinkingMatrix(cfg, grid);
            end
        end

    end
end