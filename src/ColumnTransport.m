classdef ColumnTransport
    % COLUMNTRANSPORT  1-D transport: upwind advection + flux-form diffusion.
    %
    % One explicit Euler step for a state matrix Y (n_z x n_sec).
    % Each column of Y is a depth profile for one size bin.
    %
    % Advection (sinking, downward positive z):
    %   Upwind flux at face k+1/2 = w(k) * Y(k,:)
    %   dY(k,:)/dt = -(adv_flux(k) - adv_flux(k-1)) / dz
    %   Top boundary (k=1): no inflow from above, flux(0) = 0.
    %   Bottom boundary: particles exit at k=n_z (open base).
    %
    % Diffusion (flux-form, zero-flux at both boundaries):
    %   diff_flux(k+1/2) = Kz_face(k+1/2) * (Y(k+1,:) - Y(k,:)) / dz
    %   dY(k,:)/dt = (diff_flux(k+1/2) - diff_flux(k-1/2)) / dz
    %   Kz at faces = arithmetic mean of adjacent cell values.

    methods (Static)
        function Y_new = step(Y, w_z, Kz_z, dz, dt)
            % STEP  One explicit Euler transport step.
            %
            % Inputs:
            %   Y     - n_z x n_sec, current concentrations
            %   w_z   - n_z x n_sec, sinking speed per depth & size [m/day]
            %   Kz_z  - n_z x 1, vertical diffusivity [m^2/s]
            %   dz    - scalar, cell thickness [m]
            %   dt    - scalar, time step [day]
            %
            % Output:
            %   Y_new - n_z x n_sec

            [n_z, n_sec] = size(Y);

            % ---- advection (upwind, downward positive) ----
            % face flux between layer k and k+1: use upwind (take from layer k)
            % adv_flux(k,:) = w(k,:) * Y(k,:)  for k = 1..n_z (downward face of cell k)
            adv_flux         = w_z .* max(Y, 0);  % n_z x n_sec; clip to avoid neg flux

            % net advective tendency: inflow from above minus outflow below
            % top: no inflow (open surface, particles don't come in from above)
            % bottom: open (particles leave)
            adv_tend         = zeros(n_z, n_sec);
            adv_tend(1,:)    = -adv_flux(1,:);                          % top cell: only outflow
            adv_tend(2:end,:)=  adv_flux(1:end-1,:) - adv_flux(2:end,:); % inflow - outflow

            % ---- diffusion (flux-form, [m^2/s] -> convert to [m^2/day]) ----
            day_to_sec = 8.64e4;
            Kz_day     = Kz_z * day_to_sec;  % m^2/day

            if n_z == 1
                % no interior faces when only one layer
                diff_tend = zeros(1, n_sec);
            else
                % Kz at cell faces: arithmetic mean of neighbors
                Kz_face    = 0.5 * (Kz_day(1:end-1) + Kz_day(2:end));  % (n_z-1) x 1

                % diffusive flux at interior faces
                diff_flux      = Kz_face .* (Y(2:end,:) - Y(1:end-1,:)) / dz;  % (n_z-1) x n_sec

                % add zero-flux boundaries (top and bottom face = 0)
                diff_flux_full = [zeros(1, n_sec); diff_flux; zeros(1, n_sec)];  % (n_z+1) x n_sec
                diff_tend      = (diff_flux_full(2:end,:) - diff_flux_full(1:end-1,:)) / dz;
            end

            % ---- explicit Euler update ----
            Y_new = Y + dt * (adv_tend / dz + diff_tend);

            % clip tiny negatives from numerical noise (same threshold as 0-D model)
            Y_new = max(Y_new, 0);
        end

        function cfl = maxCFL(w_z, Kz_z, dz, dt)
            % Check advective and diffusive CFL numbers.
            % Returns max across all bins and depths.
            day_to_sec = 8.64e4;
            cfl_adv  = max(abs(w_z(:))) * dt / dz;
            cfl_diff = max(Kz_z(:)) * day_to_sec * dt / dz^2;
            cfl      = max(cfl_adv, cfl_diff);
        end
    end
end
