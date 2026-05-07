classdef KernelLibrary < handle
    %KERNELLIBRARY Collection of coagulation kernel implementations
    
    methods (Static)
        function b = brownian(r, rcons, param)
            %BROWNIAN Brownian motion kernel
            % r = 2xN vector of particle radii [cm]
            % rcons and param are dummy variables for efficiency
            
            b = (2 + r(1,:)./r(2,:) + r(2,:)./r(1,:));
        end
        
        function b = curvilinearDS(r, rcons, param)
            %CURVILINEARDS Curvilinear differential sedimentation kernel
            % r = column vector of particle radii [cm]
            % rcons = column vector of conserved volume radii [cm]
            % param = parameters including r_to_rg and setcon
            
            r_small = min(r) * param.r_to_rg;
            rcons_3 = rcons .* rcons .* rcons;
            
            b = 0.5 * pi * abs(rcons_3(1,:)./r(1,:) - rcons_3(2,:)./r(2,:)) .* r_small .* r_small;
        end

        function b = curvilinearDSSinkingLaw(r, rcons, param)
            %CURVILINEARDSSINKINGLAW DS kernel using the selected sinking law.
            % Return beta/setcon so BetaAssembler can use the same scaling path.

            w = KernelLibrary.pairSettlingVelocity(r, rcons, param);
            r_small = min(r) * param.r_to_rg;
            beta_cms = 0.5 * pi * abs(w(1,:) - w(2,:)) .* r_small .* r_small;
            b = beta_cms ./ max(param.setcon, realmin);
        end
        
        function b = curvilinearShear(r, rcons, param)
            %CURVILINEARSHEAR Curvilinear shear kernel
            % r = column vector of particle radii [cm]
            % param = parameters including r_to_rg
            
            rg = (r(1,:) + r(2,:)) * param.r_to_rg;
            
            p = min(r) ./ max(r);
            p1 = 1.0 + p;
            p5 = p1 .* p1 .* p1 .* p1 .* p1;
            
            efficiency = 1.0 - (1.0 + 5.0*p + 2.5*p.*p) ./ p5;
            
            b = sqrt(8.0*pi/15.0) * efficiency .* rg .* rg .* rg;
        end
        
        function b = fractalDS(r, rcons, param)
            %FRACTALDS Fractal differential sedimentation kernel
            % r = column vector of particle radii [cm]
            % rcons = column vector of conserved volume radii [cm]
            % param = parameters including r_to_rg
            
            c1 = 0.984;  % Constant from Li and Logan
            
            rg = (r(1,:) + r(2,:)) * param.r_to_rg;
            r_ratio = min(r(1,:)./r(2,:), r(2,:)./r(1,:));
            rcons_3 = rcons .* rcons .* rcons;
            
            b = pi * abs(rcons_3(1,:)./r(1,:) - rcons_3(2,:)./r(2,:)) .* rg .* rg;
            b = b .* r_ratio.^c1;
        end
        
        function b = fractalShear(r, rcons, param)
            %FRACTALSHEAR Fractal shear kernel
            % r = column vector of particle radii [cm]
            % param = parameters including r_to_rg
            
            c1 = 0.785;  % Constant from Li and Logan
            r_ratio = min(r(1,:)./r(2,:), r(2,:)./r(1,:));
            rg = (r(1,:) + r(2,:)) * param.r_to_rg;
            
            b = 1.3 * rg .* rg .* rg;
            b = b .* r_ratio.^c1;
        end
        
        function b = rectilinearDS(r, rcons, param)
            %RECTILINEARDS Rectilinear differential sedimentation kernel
            % r = column vector of particle radii [cm]
            % rcons = column vector of conserved volume radii [cm]
            % param = parameters including r_to_rg
            
            rg = (r(1,:) + r(2,:)) * param.r_to_rg;
            rcons_3 = rcons .* rcons .* rcons;
            
            b = pi * abs(rcons_3(1,:)./r(1,:) - rcons_3(2,:)./r(2,:)) .* rg .* rg;
        end
        
        function b = rectilinearShear(r, rcons, param)
            %RECTILINEARSHEAR Rectilinear shear kernel
            % r = column vector of particle radii [cm]
            % param = parameters including r_to_rg
            
            rg = (r(1,:) + r(2,:)) * param.r_to_rg;
            b = 1.3 * rg .* rg .* rg;
        end
        
        function kernel = getKernel(kernelName)
            %GETKERNEL Get kernel function handle by name
            switch kernelName
                case 'KernelBrown'
                    kernel = @KernelLibrary.brownian;
                case 'KernelCurDS'
                    kernel = @KernelLibrary.curvilinearDS;
                case 'KernelCurDSSinkingLaw'
                    kernel = @KernelLibrary.curvilinearDSSinkingLaw;
                case 'KernelCurSh'
                    kernel = @KernelLibrary.curvilinearShear;
                case 'KernelFracDS'
                    kernel = @KernelLibrary.fractalDS;
                case 'KernelFracSh'
                    kernel = @KernelLibrary.fractalShear;
                case 'KernelRectDS'
                    kernel = @KernelLibrary.rectilinearDS;
                case 'KernelRectSh'
                    kernel = @KernelLibrary.rectilinearShear;
                otherwise
                    error('Unknown kernel: %s', kernelName);
            end
        end

        function w = pairSettlingVelocity(r, rcons, param)
            cfg = param.constants;
            law = lower(string(cfg.sinking_law));

            switch law
                case "current"
                    r_v = rcons;
                    r_i = KernelLibrary.conservativeToFractalRadius(r_v, cfg);
                    w = param.setcon .* (r_v .^ 3) ./ max(r_i, realmin);

                case "kriest_8"
                    d_cm = KernelLibrary.pairDiameterCm(r, rcons, param);
                    w = KernelLibrary.mdayToCms(66 .* d_cm .^ 0.62, cfg);

                case "kriest_9"
                    d_cm = KernelLibrary.pairDiameterCm(r, rcons, param);
                    w = KernelLibrary.mdayToCms(132 .* d_cm .^ 0.62, cfg);

                case "siegel_2025"
                    d_cm = KernelLibrary.pairDiameterCm(r, rcons, param);
                    d_mm = d_cm * 10.0;
                    w = KernelLibrary.mdayToCms(20.2 .* d_mm .^ 0.67, cfg);

                case "kriest_8_capped"
                    d_cm = KernelLibrary.pairDiameterCm(r, rcons, param);
                    w_mday = 66 .* d_cm .^ 0.62;
                    w_max = cfg.sinking_w_max_mday;
                    w = KernelLibrary.mdayToCms(min(w_mday, w_max), cfg);

                case "kriest_8_flat"
                    d_cm = KernelLibrary.pairDiameterCm(r, rcons, param);
                    w_mday = 66 .* d_cm .^ 0.62;
                    d_flat_cm = cfg.sinking_d_flat_cm;
                    w_flat = 66 .* d_flat_cm .^ 0.62;
                    w_mday(d_cm >= d_flat_cm) = w_flat;
                    w = KernelLibrary.mdayToCms(w_mday, cfg);

                otherwise
                    error("Unknown sinking_law: %s", cfg.sinking_law);
            end

            if isprop(cfg, 'sinking_scale') && ~isempty(cfg.sinking_scale)
                w = w * cfg.sinking_scale;
            end
            w(~isfinite(w)) = 0;
            w(w < 0) = 0;
        end

        function d_cm = pairDiameterCm(r, rcons, param)
            cfg = param.constants;
            if isprop(cfg, 'sinking_size') && strcmpi(cfg.sinking_size, 'image')
                d_cm = 2.0 * param.r_to_rg .* r;
            else
                d_cm = 2.0 .* rcons;
            end
        end

        function v_cms = mdayToCms(w_mday, cfg)
            v_cms = (w_mday * 100) / cfg.day_to_sec;
        end

        function amfrac = currentAmfrac(cfg)
            a0 = cfg.d0 / 2.0;
            amfrac_temp = (4.0 / 3.0 * pi) ^ (-1.0 / cfg.fr_dim) * a0 ^ (1.0 - 3.0 / cfg.fr_dim);
            amfrac = amfrac_temp * sqrt(0.6);
        end

        function setcon = currentSetcon(cfg)
            del_rho = (4.5 * 2.48) * cfg.kvisc * cfg.rho_fl / cfg.g * (cfg.d0 / 2.0) ^ (-0.83);
            setcon = (2.0 / 9.0) * del_rho / cfg.rho_fl * cfg.g / cfg.kvisc;
        end

        function r_i = conservativeToFractalRadius(r_v, cfg)
            amfrac = KernelLibrary.currentAmfrac(cfg);
            bmfrac = 1.0 / cfg.fr_dim;
            av_vol = (4.0 / 3.0) .* pi .* (r_v .^ 3);
            r_i = amfrac .* (av_vol .^ bmfrac);
        end

        function r_v = fractalToConservativeRadius(r_i, cfg)
            amfrac = KernelLibrary.currentAmfrac(cfg);
            av_vol = (r_i ./ max(amfrac, realmin)) .^ cfg.fr_dim;
            r_v = ((0.75 / pi) .* av_vol) .^ (1.0 / 3.0);
        end
    end
end
