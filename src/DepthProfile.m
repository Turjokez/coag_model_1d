classdef DepthProfile
    % DEPTHPROFILE  Ocean physics profiles on a depth grid.
    %
    % Stores T, S, rho, kinematic viscosity, turbulent dissipation,
    % and diffusivity at each depth cell center.
    % All profiles are in consistent units (see property comments).
    %
    % Usage:
    %   p = DepthProfile.typical(z_centers);  % simple default profile
    %   p = DepthProfile(z_centers, T_K, nu, eps, Kz);

    properties
        z       % depth cell centers [m], n_z x 1
        T_K     % temperature [K],             n_z x 1
        S       % salinity [psu],               n_z x 1
        rho     % density [g/cm^3],             n_z x 1
        nu      % kinematic viscosity [cm^2/s], n_z x 1
        eps     % turbulent dissipation [cm^2/s^3], n_z x 1
        Kz      % vertical diffusivity [m^2/s], n_z x 1
        Zc      % filter feeder concentration [m^-3], n_z x 1
        Zf      % flux feeder concentration [m^-3], n_z x 1
    end

    methods
        function obj = DepthProfile(z, T_K, S, rho, nu, eps, Kz)
            % Constructor: pass all profiles directly.
            obj.z   = z(:);
            obj.T_K = T_K(:);
            obj.S   = S(:);
            obj.rho = rho(:);
            obj.nu  = nu(:);
            obj.eps = eps(:);
            obj.Kz  = Kz(:);
            obj.Zc  = [];
            obj.Zf  = [];
        end

        function scale = brownianScale(obj, cfg)
            % Ratio of Brownian prefactor at each depth to the reference value.
            % Reference is grid.conBr = 2*k_B*T_ref / (3*mu_ref).
            % At depth k: conBr(k) = 2*k_B*T_K(k) / (3*mu(k))
            %             mu(k) = rho(k) * nu(k)
            % scale(k) = conBr(k) / conBr_ref
            mu_ref  = cfg.rho_fl * cfg.kvisc;    % g/(cm*s)
            mu_z    = obj.rho .* obj.nu;          % g/(cm*s)
            scale   = (obj.T_K / cfg.temp) .* (mu_ref ./ mu_z);
        end

        function scale = shearScale(obj, cfg)
            % Ratio of shear rate G(k) = sqrt(eps(k)/nu(k)) to reference G_ref.
            % Reference G_ref = cfg.gamma [s^-1].
            % eps in cm^2/s^3, nu in cm^2/s -> G in s^-1.
            G_z    = sqrt(obj.eps ./ obj.nu);    % s^-1 per depth
            G_ref  = cfg.gamma;                  % s^-1 (reference)
            scale  = G_z / max(G_ref, 1e-30);
        end

        function scale = dsScale(obj, cfg)
            % Ratio of sinking speed at each depth to reference.
            % w(k) = w_ref * nu_ref / nu(k)  (viscosity correction)
            % DS kernel scales linearly with sinking speed.
            scale = cfg.kvisc ./ obj.nu;
        end
    end

    methods (Static)
        function p = typical(z_centers)
            % Simple open-ocean profile with Stemmann zoo profiles.
            z = z_centers(:);
            H = max(z);
            if H < 1, H = 2000; end

            T_C = 20 - 16 * z / H;
            T_K = T_C + 273.15;
            S   = 35 * ones(size(z));
            rho = 1.025 + 0.0005 * z / H;
            nu  = 0.01 * exp(0.02 * (20 - T_C));

            eps_surf = 1e-4;
            H_mix    = 100;
            eps = eps_surf * exp(-z / H_mix) + 1e-8;

            Kz = 1e-5 * ones(size(z));
            Kz(z <= H_mix) = 1e-3;

            p = DepthProfile(z, T_K, S, rho, nu, eps, Kz);

            % add depth-varying zoo from Stemmann 2004 Fig 1
            [p.Zc, p.Zf] = DepthProfile.stemmannZoo(z);
        end

        function [Zc, Zf] = stemmannZoo(z)
            % Stemmann 2004 Fig 1 zooplankton profiles.
            % Relative abundance x max concentration for each group.
            % Filter feeders:  max = 0.307 m^-3, peak ~350 m
            % Flux feeders:    max = 0.063 m^-3, increases with depth
            z = z(:);

            % control points from Fig 1 [depth_m, relative_abundance]
            z_pts = [0; 100; 200; 300; 400; 500; 600; 700; 800; 900; 1000];

            ff_rel  = [0.10; 0.30; 0.70; 0.95; 1.00; 0.85; 0.60; 0.40; 0.25; 0.15; 0.10];
            flx_rel = [0.05; 0.10; 0.20; 0.30; 0.40; 0.50; 0.60; 0.70; 0.80; 0.90; 1.00];

            % interpolate to model depths, clamp to [0,1]
            ff  = interp1(z_pts, ff_rel,  z, 'pchip', 'extrap');
            flx = interp1(z_pts, flx_rel, z, 'pchip', 'extrap');
            ff  = max(0, min(1, ff));
            flx = max(0, min(1, flx));

            Zc = 0.307 * ff;   % filter feeders [m^-3]
            Zf = 0.063 * flx;  % flux feeders   [m^-3]
        end
    end
end
