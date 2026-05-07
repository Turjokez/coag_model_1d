classdef ColumnGrid
    % COLUMNGRID  Uniform depth grid for the 1-D column model.
    %
    % Divides [0, H] into n_z equal cells of thickness dz = H/n_z.
    % Cell centers are at dz/2, 3*dz/2, ..., H - dz/2.
    % Cell faces are at 0, dz, 2*dz, ..., H.

    properties
        n_z       % number of depth cells
        H         % total column depth [m]
        dz        % cell thickness [m]
        z_edges   % (n_z+1) x 1 face depths [m]
        z_centers % n_z x 1 cell center depths [m]
    end

    methods
        function obj = ColumnGrid(H, n_z)
            % ColumnGrid(H, n_z): H in meters, n_z integer.
            obj.H         = H;
            obj.n_z       = n_z;
            obj.dz        = H / n_z;
            obj.z_edges   = (0 : n_z)' * obj.dz;
            obj.z_centers = obj.z_edges(1:end-1) + obj.dz / 2;
        end
    end
end
