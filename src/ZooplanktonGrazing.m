%> @brief Zooplankton grazing module (Stemmann et al. 2004).
%> @details Implements filter-feeder clearance and flux-feeder interception
%>          with fecal pellet production at a specified size bin.
classdef ZooplanktonGrazing
    %ZOOPLANKTONGRAZING Stemmann-style zooplankton grazing for 0-D slab.

    properties
        % Filter feeders
        Zc = 100;      % filter feeder concentration [ind m^-3]
        c  = 1e-4;     % clearance rate per individual [m^3 ind^-1 day^-1]

        % Flux feeders
        Zf = 50;       % flux feeder concentration [ind m^-3]
        s  = 1e-4;     % cross-sectional capture area [m^2 ind^-1]

        % Fecal production
        p  = 0.3;      % egestion fraction
        ic = 1;        % minimum fecal pellet section index
    end

    methods
        function obj = ZooplanktonGrazing(varargin)
            % Constructor: name-value pairs override defaults.
            if nargin > 0
                for k = 1:2:numel(varargin)
                    if isprop(obj, varargin{k})
                        obj.(varargin{k}) = varargin{k+1};
                    end
                end
            end
        end

        %> @brief Compute grazing loss rate and fecal pellet production.
        %> @param v      Biovolume spectrum [cm^3 cm^-3], length n_sec.
        %> @param w_cms  Sinking speed array [cm s^-1], length n_sec.
        %> @param Zc_in  Filter feeder concentration override [ind m^-3] (optional).
        %> @param Zf_in  Flux feeder concentration override [ind m^-3] (optional).
        %> @return dvdt    Loss rate [bv day^-1], length n_sec.
        %> @return fp_flux Fecal pellet production [bv day^-1], length n_sec.
        function [dvdt, fp_flux] = graze(obj, v, w_cms, Zc_in, Zf_in)
            % GRAZE  Grazing tendency following Stemmann et al. 2004.
            % Optional Zc_in, Zf_in override obj.Zc, obj.Zf (for depth profiles).
            %
            % Inputs:
            %   v      - n_sec x 1 biovolume concentration
            %   w_cms  - n_sec x 1 settling velocity [cm/s]
            % Outputs:
            %   dvdt    - n_sec x 1 tendency [per day] applied to aggregate array
            %   fp_flux - scalar [bv day^-1], total fecal production at this layer
            %
            % When called with one output (old usage), fecal return is added back
            % into dvdt at bin ic+1 (backward compatible).
            % When called with two outputs, fecal is returned separately so the
            % caller can route it to a separate fecal pellet array (Y_fp).
            v     = v(:);
            w_cms = w_cms(:);
            n     = numel(v);

            % use depth-specific values if provided
            if nargin >= 4 && ~isempty(Zc_in), Zc = Zc_in; else, Zc = obj.Zc; end
            if nargin >= 5 && ~isempty(Zf_in), Zf = Zf_in; else, Zf = obj.Zf; end

            % Convert settling speed to m/day.
            day_to_sec = 8.64e4;
            w_mday = (w_cms / 100) * day_to_sec;

            % Removal rates [day^-1].
            rate_FF = obj.c * Zc;
            rate_FL = w_mday * obj.s * Zf;
            rate    = rate_FF + rate_FL;

            % Consumption per bin.
            consumption = rate .* v;
            actual_consumption = min(consumption, v);

            % dvdt: only losses from the aggregate array.
            dvdt = -actual_consumption;

            % Total fecal production [bv day^-1].
            fp_flux = obj.p * sum(actual_consumption);

            if nargout < 2
                % Old usage: put fecal return back into aggregate array (bin ic+1).
                target_bin = max(1, min(n, round(obj.ic) + 1));
                if fp_flux > 0
                    dvdt(target_bin) = dvdt(target_bin) + fp_flux;
                end
            end
            % New usage (nargout == 2): caller gets fp_flux separately and
            % routes it to Y_fp. Nothing added to dvdt here.
        end

        function [dvdt, fp_flux] = mine(obj, v, w_cms, av_vol, Zm, dm_gut, s_area, min_bin)
            % MINE  Partial consumption (mining) following Stemmann 2004 Part I Eq. 25.
            %
            % Small copepods take a fixed bite dm_gut from each particle contact.
            % Particles shrink from bin i toward bin i-1 (not removed entirely).
            % Mining only applied to bins >= min_bin (particles large enough to mine).
            % Fecal fraction p of consumed mass goes to fp_flux.
            %
            % Inputs:
            %   v        - n_sec x 1 biovolume concentration
            %   w_cms    - n_sec x 1 settling velocity [cm/s]
            %   av_vol   - n_sec x 1 average particle biovolume per bin [cm^3]
            %   Zm       - miner concentration [ind m^-3]
            %   dm_gut   - mass per contact [cm^3] (gut volume)
            %   s_area   - cross-section area [m^2 ind^-1]
            %   min_bin  - first bin where mining is active (default 12, ~254 um)
            %
            % Outputs:
            %   dvdt    - n_sec x 1 tendency [bv day^-1] for aggregate array
            %   fp_flux - scalar [bv day^-1], fecal production from mining
            if nargin < 8 || isempty(min_bin), min_bin = 12; end

            v      = v(:);
            w_cms  = w_cms(:);
            av_vol = av_vol(:);
            n      = numel(v);

            day_to_sec = 8.64e4;
            w_mday = (w_cms / 100) * day_to_sec;  % m/day

            % encounter rate per bin [day^-1]: same flux-feeding structure
            em = w_mday * s_area * Zm;

            % bite rate: dm_gut limited to particle volume
            dm_eff   = min(dm_gut, av_vol);
            bite_fac = em .* dm_eff ./ av_vol;

            % zero out bins below min_bin — small particles are not mined
            bite_fac(1:min_bin-1) = 0;

            dvdt = zeros(n, 1);
            for i = 1:n
                % two losses from bin i: bite mass + shrunken particle exits bin
                dvdt(i) = dvdt(i) - 2 * bite_fac(i) * v(i);
                % gain from bin i+1 shrinking in
                if i < n
                    dvdt(i) = dvdt(i) + bite_fac(i+1) * v(i+1);
                end
            end

            % fecal: fraction p of total consumed mass (active bins only)
            fp_flux = obj.p * sum(bite_fac .* v);
        end
    end
end
