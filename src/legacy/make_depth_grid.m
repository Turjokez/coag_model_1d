function grid = make_depth_grid(z_max_m, dz_m)
% make_depth_grid
% Build simple depth grid structure.

grid = struct();
grid.z_max_m = z_max_m;
grid.dz_m = dz_m;
grid.z_m = (0:dz_m:z_max_m)';
grid.nz = numel(grid.z_m);
end

