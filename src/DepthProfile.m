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
            % TYPICAL  Simple but realistic open-ocean profile.
            % T: linear 20 C at surface to 4 C at 2000 m
            % S: constant 35 psu
            % rho: 1.025 + small depth gradient
            % nu: increases with depth as water gets colder
            % eps: exponential decay with scale depth 100 m
            % Kz: large in mixed layer (top 100 m), small below
            z = z_centers(:);
            H = max(z);
            if H < 1, H = 2000; end

            % temperature: linear from 20 C at surface to 4 C at bottom
            T_C = 20 - 16 * z / H;
            T_K = T_C + 273.15;

            % salinity: constant
            S = 35 * ones(size(z));

            % density [g/cm^3]: simple linear approximation
            rho = 1.025 + 0.0005 * z / H;

            % kinematic viscosity [cm^2/s]:
            % nu ~ 0.01 at 20 C, increases ~37% at 4 C
            % simple exponential fit: nu = 0.01 * exp(0.02*(20-T_C))
            nu = 0.01 * exp(0.02 * (20 - T_C));

            % turbulent dissipation [cm^2/s^3]:
            % surface: 1e-4 cm^2/s^3 (= 1e-8 W/kg), decays with 100 m scale
            % 1 W/kg = 1 m^2/s^3 = 1e4 cm^2/s^3
            eps_surf     = 1e-4;      % cm^2/s^3
            H_mix        = 100;       % m
            eps = eps_surf * exp(-z / H_mix) + 1e-8;  % floor at 1e-8

            % vertical diffusivity [m^2/s]:
            % mixed layer (~100 m): 1e-3 m^2/s, below: 1e-5 m^2/s
            Kz = 1e-5 * ones(size(z));
            Kz(z <= H_mix) = 1e-3;

            p = DepthProfile(z, T_K, S, rho, nu, eps, Kz);
        end
    end
end
