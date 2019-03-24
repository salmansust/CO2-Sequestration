sector        = 'NPD5';
grdecl        = readGRDECL([sector, '.grdecl']);
actnum        = grdecl.ACTNUM;
grdecl.ACTNUM = ones(prod(grdecl.cartDims),1);
G             = processGRDECL(grdecl, 'checkgrid', false);

grdecl.ACTNUM = actnum; clear actnum;
G = processGRDECL(grdecl); clear grdecl;
G = computeGeometry(G);

% Plotting a height map of the field using the z-component of the centroids
% of the cells
clf,
plotCellData(G,G.cells.centroids(:,3),'EdgeColor','k','EdgeAlpha',0.1);
colorbar, view(3), axis tight off, view(-20,40), zoom(1.2)

