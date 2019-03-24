% The following line ensures that we have the necessary modules loaded
moduleCheck co2lab

% Loading the grid from disk
coarsening_level = 1; % grid downsampling factor (1 = no downsampling)
Gt = getFormationTopGrid('Johansenfm', coarsening_level);

% We check the number of cells in the top surface grid.
Gt.cells.num

% We then perform the trapping analysis.
ta = trapAnalysis(Gt, false)

% Top view
figure
plotCellData(Gt, ta.traps, 'edgealpha', 0.1);
view(0, 90); axis tight; colormap lines;
set(gcf, 'position', [10 10 500 800]);

% Oblique view
figure
plotCellData(Gt, ta.traps, 'edgealpha', 0.1);
view(290, 60); axis tight; colormap lines;
set(gcf, 'position', [520 10 1200 800]);

trapfield = ta.traps;
trapfield(trapfield==0) = NaN;

% Top view
figure;
plotCellData(Gt, trapfield, 'edgecolor', 'none');
plotCellData(Gt, ta.trap_regions, 'facealpha', 0.3, 'edgecolor', 'none');
view(0, 90); axis tight; colormap lines;
set(gcf, 'position', [10 10 500 800]);

% Oblique view
figure;
plotCellData(Gt, trapfield, 'edgecolor', 'none');
plotCellData(Gt, ta.trap_regions, 'facealpha', 0.3, 'edgealpha', 0.1);
view(290, 60); axis tight; colormap lines;
set(gcf, 'position', [520 10 1200 800]);

% We compute the number of cells in each trap region, and store the number in
% 'region_cellcount'.  This variable is thus a vector where entry 'i'
% states the number of cells belonging to trap region 'i'.
% We then produce a bar plot from the sorted result.
region_cellcount = ...
    diff(find(diff([sort(ta.trap_regions); max(ta.trap_regions)+1])));
bar(sort(region_cellcount, 'descend'));
hold on;

% We then do the same for traps.
trap_cellcount = diff(find(diff([sort(ta.traps); max(ta.traps)+1])));
bar(sort(trap_cellcount, 'descend'), 'r');

% Find index of trap with the largest spill region
[~, ix] = max(region_cellcount);

% Visualize the trap along with the spill region
figure;
plotGrid(topSurfaceGrid(extractSubgrid(Gt.parent, find(ta.traps==ix))), ...
         'facecolor', 'r', ...
         'edgealpha', 0.1)
plotGrid(topSurfaceGrid(extractSubgrid(Gt.parent, find(ta.trap_regions==ix))), ...
         'facecolor','r', ...
         'facealpha', 0.1, ...
         'edgealpha', 0.1);
view(-64, 36);
set(gcf, 'position', [10 10 1000 700]);

% Identifying "river" cells
river_field = NaN(size(ta.traps));
for i = 1:numel(ta.cell_lines)
   rivers = ta.cell_lines{i};
   for r = rivers
      river_field(r{:}) = 1;
   end
end

% Producing figure with traps and rivers (top view)
figure;
plotCellData(Gt, river_field, 'edgealpha', 0.1);
plotCellData(Gt, trapfield, 'edgecolor', 'none');
view(0, 90); axis tight; colormap lines;
set(gcf, 'position', [10 10 500 800]);

% Same figure, but oblique view
figure;
plotCellData(Gt, river_field, 'edgealpha', 0.1);
plotCellData(Gt, trapfield, 'edgecolor', 'none');
view(290, 60); axis tight; colormap lines;
set(gcf, 'position', [520 10 1200 800]);

% Producing a topographical map of caprock surface, traps and rivers
h = figure;
mapPlot(h, Gt, 'traps', ta.traps, 'rivers', ta.cell_lines);
set(gcf, 'position', [10 10 500 800]);

% We add a zero to the vector of spill point depths.  This is necessary for
% indexing purposes, as explained in the next comment.
depths = [0; ta.trap_z];

% The line below gives the correct "trap height" value for cells within
% traps.  For the other cells, the value is incorrect, but will be fixed by
% the succeeding line.  Note that since we have added an entry to the front
% of the 'depths' vector, we increase indices by 1.  (Otherwise, the
% presence of zeros in 'ta.traps' would cause an indexing error, since Matlab
% arrays are indexed from 1, not 0).
h = min((depths(ta.traps+1) - Gt.cells.z), Gt.cells.H);
h(ta.traps==0) = 0; % cells outside any trap should have zero trap volume

% Bulk cell volumes are now simply obtained by multiplying by the respective
% areas.  For non-trap cells, the volumes will be zero, since 'h' is zero for
% these cells.
cell_tvols = h .* Gt.cells.volumes;

% We accumulate the trap volume of individual cells to obtain the total bulk
% trap volume for each structural trap.
tvols = accumarray(ta.traps+1, cell_tvols);

% The first entry of 'tvols' contain the combined trap volume of all non-trap
% cells.  Obviously, this value should be zero, but let us check this to make
% sure.
assert(tvols(1)==0);

% We remove the first entry of 'tvols', since we are only interested in the
% volumes of the actual traps.
tvols = tvols(2:end);

% We produce a bar plot where trap volumes are plotted in descending order.
% Here too, we notice the presence of a handful of large traps, and a long
% tail of vanishingly small traps.
figure;
bar(sort(tvols, 'descend'))

porosity = 0.1071; % (This porosity value for Statfjord is from the Norwegian
                   % Petroleum Directorate)
seafloor_temp = 7 + 273.15; % Seafloor temperature, in degrees Kelvin
temp_grad = 30; % temperature increase per kilometer depth
rho_brine = 1000; % brine density (kilogram per cubic meter)

% Ensure that gravity is not zero
gravity on;

% computing temperature field (at caprock level)
T = seafloor_temp + temp_grad .* Gt.cells.z/1000;

% computing hydrostatic pressure (at caprock level)
P = 1 * atm + rho_brine * norm(gravity) * Gt.cells.z;

% Making stripped-down fluid object that only contains CO? density function
fluid = addSampledFluidProperties(struct, 'G');

% Computing local CO? (at caprock level)
CO2_density = fluid.rhoG(P, T);

% Computing trap capacity in mass term for each cell
cell_tmass = cell_tvols .* porosity .* CO2_density;

% Accumulating cell values to get trap capacity in mass term for each trap
tmass = accumarray(ta.traps+1, cell_tmass);
assert(tmass(1)==0);
% Normal view
figure
plotCellData(Gt, tmass(ta.traps+1), 'edgealpha', 0.1);
view(0, 90); axis tight; colormap cool; colorbar;
set(gcf, 'position', [10 10 500 800]);

% Oblique view
figure
plotCellData(Gt, tmass(ta.traps+1), 'edgealpha', 0.1);
view(290, 60); axis tight; colormap cool; colorbar;
set(gcf, 'position', [520 10 1200 800]);

% We remove the first entry of 'tmass' to obtain a vector where each entry
% represents the structural capacity in mass terms of the trap with
% correspodning index. (The removed element represents the region outside any
% spill region, so it has a trapping capacity of zero).
tmass = tmass(2:end);

cum_reachable = zeros(size(ta.traps));

% Adding the capacity of each trap to its own spill region and all downstream
% regions
for trap_ix = 1:max(ta.traps)

   region = ta.trap_regions==trap_ix;

   % Counting this trap's capacity towards all cells in its spill region
   cum_reachable(region) = cum_reachable(region) + tmass(trap_ix);

   visited_regions = trap_ix;

   % Computing contribution to cells associated with downstream traps
   downstream = find(ta.trap_adj(:,trap_ix));
   while ~isempty(downstream)

      region = ta.trap_regions == downstream(1);
      cum_reachable(region) = cum_reachable(region) + tmass(trap_ix);

      visited_regions = [visited_regions;downstream(1)];

      % downstream(1) has now been processed, so we remove it from the
      % downstream vector.
      downstream = [downstream(2:end); find(ta.trap_adj(:,downstream(1)))];
      downstream = setdiff(downstream, visited_regions);
   end
end

% we divide 'cum_reachable' by 1e12 to have the result in gigatons
cum_reachable = cum_reachable / 1e12;

% Plot result (normal view)
figure
plotCellData(Gt, cum_reachable, 'edgealpha', 0.1);
view(0, 90); axis tight; colormap cool; colorbar;
set(gcf, 'position', [10 10 500 800]);

% Plot result (oblique view)
figure
plotCellData(Gt, cum_reachable, 'edgealpha', 0.1);
view(290, 60); axis tight; colormap cool; colorbar;
set(gcf, 'position', [520 10 1200 800]);
