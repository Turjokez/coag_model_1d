classdef SimulationConfig < matlab.mixin.Copyable
    %SIMULATIONCONFIG Configuration class for coagulation simulation parameters
    
    properties
        % Physical parameters
        rho_fl = 1.0275;        % Fluid density [g cm^{-3}]
        kvisc = 0.01;           % Kinematic viscosity [cm^2 s^{-1}]
        g = 980;                % Accel. due to gravity [cm s^{-2}]
        day_to_sec = 8.64e04;   % Seconds in a day [s d^{-1}]
        k = 1.3e-16;           % Boltzmanns constant [erg K^{-1}]
        r_to_rg = 1.6;         % Interaction radius to radius of gyration
        box_depth = [];   % meters. If set, use w/box_depth for 0-D sinking loss
        enable_pp = false;      % on/off
        pp_bin    = 1;          % which bin gets the source
        pp_source = 0;          % source strength (state units per day)

        % Section/coagulation related parameters
        n_sections = 20;        % Number of sections
        kernel = 'KernelBrown'; % Kernel type
        d0 = 20e-4;            % Diameter of unit particle [cm]
        fr_dim = 2.33;         % Particle fractal dimension
        n1 = 100;              % No. particles cm^{-3} in first section
        
        sinking_law  = 'current';
        sinking_size = 'volume';
        ds_kernel_mode = 'sinking_law'; % 'legacy' or 'sinking_law'
        enable_coag    = true;
        enable_sinking = true;
        enable_disagg  = false;
        enable_linear  = false;
        sinking_scale = 1.0;   % multiplier for settling velocity (dimensionless)
        sinking_w_max_mday = 70; % optional cap for capped sinking laws [m/day]
        sinking_d_flat_cm = 0.1; % optional flat threshold for kriest_8_flat [cm]

        % Other input parameters
        temp = 20 + 273;       % Temperature [K]
        alpha = 1.0;           % Stickiness
        dz = 65;               % Layer thickness [m]
        gamma = 0.1;           % Average shear rate [s^{-1}]
        growth = 0.15;         % Specific growth rate in first section [d^{-1}]
        gro_sec = 4;           % Section at which growth in aggregates starts
        num_1 = 10^3;          % Number of particle cm^{-3} in first section
        
        % Kernel component scaling (for diagnostics)
        scale_brown = 1.0;     % Brownian kernel multiplier
        scale_shear = 1.0;     % Shear kernel multiplier
        scale_ds    = 1.0;     % Differential settling kernel multiplier
        
        % Disaggregation parameters
        c3 = 0.02;             % For curvilinear kernel (disagg strength)
        c4 = 1.45;             % For curvilinear kernel
        disagg_mode = 'legacy';     % 'legacy' or 'operator_split'
        disagg_outer_dt = 1/24;     % days (1 hour) for operator-split loop
        disagg_frac_next = 2/3;     % fraction to next smaller bin
        disagg_C = 3.0;             % D_max = C * epsilon^(-gamma), C in mm
        disagg_gamma = 0.15;        % exponent for D_max scaling
        disagg_epsilon = [];        % turbulence dissipation rate (scalar or function handle)
        disagg_dmax_cm = [];        % optional override for D_max in cm
        disagg_dmax_A  = [];        % Parker calibration constant [m]; default 9.39e-6
        disagg_dmax_cap_cm = Inf;   % hard upper cap on D_max [cm]; Inf = no cap (default)
        
        % Parameters for solving equations
        t_init = 0.0;          % Initial time for integrations [d]
        t_final = 30.0;        % Final time for integrations [d]
        delta_t = 1.0;         % Time interval for output [d]
        
        % Code Runtime Options
        tracer = false;         % Integrate tracer as well [false=no, true=yes]
        proc_substeps = 10;     % substeps for explicit process-rate update in ColumnRHS

        % Zooplankton grazing (Stemmann 2004)
        enable_zoo  = false;   % turn grazing on/off
        zoo_Zc      = 0.307;   % filter feeder max [m^-3], Stemmann 2004
        zoo_c       = 1e-4;    % clearance rate [m^3 ind^-1 day^-1]
        zoo_Zf      = 0.063;   % flux feeder max [m^-3], Stemmann 2004
        zoo_s       = 1e-4;    % capture cross-section [m^2 ind^-1]
        zoo_p       = 0.3;     % egestion fraction
        zoo_ic      = 7;       % fecal pellet target = bin 8 (~115 um)

        % Fecal pellet sinking (Stokes law with dense pellet excess density)
        fp_excess_density = 0.15;  % excess density [g/cm^3] (~150 kg/m^3, cylindrical pellets)

        % Fecal pellet cross-coagulation stickiness with marine snow
        % Fecal pellets are compact (fractal dim ~2.8), less TEP -> less sticky.
        % No direct measurement exists; 0.5 is a conservative starting value.
        fp_alpha_cross = 0.5;

        % Microbial remineralization (first-order loss on Y and Y_fp)
        % Applied as operator-split exact decay: Y *= exp(-r*dt)
        enable_microbe     = false;  % turn on/off
        microbe_r0         = 0.03;   % base rate [day^-1] (Iversen & Ploug 2013)
        microbe_fp_mult    = 1.0;    % r multiplier for fecal vs aggregate (e.g. 1.5)
        microbe_use_temp   = false;  % Q10 temperature scaling on/off
        microbe_q10        = 2.0;    % Q10 factor (Iversen & Ploug 2013)
        microbe_tref_C     = 20;     % reference temperature [deg C]
        microbe_gamma_size = 0.0;    % size exponent: r ~ d^(-gamma), 0 = size-independent
        microbe_dref_cm    = 0.01;   % reference diameter for size scaling [cm] (= 100 um)

        % Micro-zooplankton mining (Stemmann 2004 Part I, Eq. 25)
        % Small copepods (e.g. Oncaea) bite a fixed chunk dm from each particle.
        % Particles shrink bin-by-bin; fecal fraction p goes to Y_fp.
        enable_mining    = false;   % on/off
        mining_Zm        = 250;     % miner concentration [ind m^-3] (Stemmann Fig 1)
        mining_dm        = 1e-5;    % mass uptake per contact [cm^3] (gut volume)
        mining_s         = 1.3e-5;  % cross-section area [m^2 ind^-1]
        mining_min_bin   = 12;      % only mine bins >= this (~254 um, marine snow size)

        % Surface production (1-D only — applies to layer 1, bin surface_pp_bin)
        enable_surface_pp   = false;   % surface phytoplankton source on/off
        surface_pp_bin      = 1;       % which size bin gets the source
        surface_pp_rate     = 1e-8;    % constant source strength [bv day^-1]
        surface_pp_mu       = 0;       % growth rate [day^-1]; if > 0, uses mu*phi mode instead

        % DVM fecal rerouting -- commented out June 16 2026
        % All test variants returned null (no change in deep ratio at 475 m).
        % Archibald 2019 rerouting, gut-memory (G_gut), and Bianchi-lite (G_gut + B_mig)
        % all tested. Params preserved for reference if Adrian asks to revisit.
        % enable_dvm      = false;
        % dvm_p           = 0.5;    % fraction of grazers that migrate [0-1]
        % dvm_ffec        = 0.7;    % fraction of migrator fecal in feeding zone [0-1]
        % dvm_feed_zmax   = 150;    % bottom of feeding zone [m]
        % dvm_zmin        = 300;    % top of deep injection band [m]
        % dvm_zmax        = 500;    % bottom of deep injection band [m]
        % dvm_tau_gut     = 0.25;   % gut clearance time [day]
        % enable_dvm_bianchi = false;
        % dvm_tau_body       = 1.0;
        % dvm_body_frac      = 0.3;
        % dvm_mort_rate      = 0.02;
        % dvm_mort_bin       = 15;
    end
    
    methods
        function obj = SimulationConfig(varargin)
            % SimulationConfig: Constructor - allows for initialization with parameter-value pairs.
            % This method enables the user to create an instance of the
            % SimulationConfig class and, optionally, override the default
            % property values by passing in a list of parameter-value pairs.
            if nargin > 0
                for i = 1:2:length(varargin)
                    if isprop(obj, varargin{i})
                        obj.(varargin{i}) = varargin{i+1};
                    end
                end
            end
        end
        
        function grid = derive(obj)
            % derive: Computes and returns a DerivedGrid object with precomputed values.
            % This method acts as a factory, creating a new object that
            % contains values derived from the simulation configuration.
            grid = DerivedGrid(obj);
        end
        
        function validate(obj)
            % validate: Validates configuration parameters to ensure they are valid for the simulation.
            % This method checks for common errors in the configuration,
            % such as non-positive values for parameters that must be
            % positive, and asserts that the specified values are logical.
            assert(obj.n_sections > 0, 'n_sections must be positive');
            assert(obj.t_final > obj.t_init, 't_final must be greater than t_init');
            assert(obj.delta_t > 0, 'delta_t must be positive');
            if isprop(obj,'sinking_law') && ~isempty(obj.sinking_law)
                law = lower(string(obj.sinking_law));
                valid_laws = ["current","kriest_8","kriest_9","siegel_2025","kriest_8_capped","kriest_8_flat"];
                assert(any(law == valid_laws), 'sinking_law must be one of: current, kriest_8, kriest_9, siegel_2025, kriest_8_capped, kriest_8_flat');
            end
            if isprop(obj,'sinking_size') && ~isempty(obj.sinking_size)
                sz = lower(string(obj.sinking_size));
                valid_sizes = ["volume","image"];
                assert(any(sz == valid_sizes), 'sinking_size must be: volume or image');
            end
            if isprop(obj,'ds_kernel_mode') && ~isempty(obj.ds_kernel_mode)
                mode = lower(string(obj.ds_kernel_mode));
                valid_modes = ["legacy","sinking_law"];
                assert(any(mode == valid_modes), 'ds_kernel_mode must be legacy or sinking_law');
            end
            if isprop(obj,'sinking_scale') && ~isempty(obj.sinking_scale)
                assert(isfinite(obj.sinking_scale) && obj.sinking_scale >= 0, 'sinking_scale must be finite and >= 0');
            end
            if isprop(obj,'disagg_mode') && ~isempty(obj.disagg_mode)
                mode = lower(string(obj.disagg_mode));
                valid_modes = ["legacy","operator_split","logistic"];
                assert(any(mode == valid_modes), 'disagg_mode must be legacy, operator_split, or logistic');
            end
            if isprop(obj,'disagg_outer_dt') && ~isempty(obj.disagg_outer_dt)
                assert(isfinite(obj.disagg_outer_dt) && obj.disagg_outer_dt > 0, 'disagg_outer_dt must be > 0');
            end
            if isprop(obj,'disagg_frac_next') && ~isempty(obj.disagg_frac_next)
                assert(obj.disagg_frac_next >= 0 && obj.disagg_frac_next <= 1, 'disagg_frac_next must be in [0,1]');
            end
            if isprop(obj,'enable_disagg') && obj.enable_disagg
                if isprop(obj,'disagg_mode') && strcmpi(string(obj.disagg_mode), 'operator_split')
                    has_dmax = isprop(obj,'disagg_dmax_cm') && ~isempty(obj.disagg_dmax_cm);
                    has_eps  = isprop(obj,'disagg_epsilon') && ~isempty(obj.disagg_epsilon);
                    assert(has_dmax || has_eps, 'operator_split disagg needs disagg_dmax_cm or disagg_epsilon');
                end
            end
            % Add more validation as needed
        end

    end
end
