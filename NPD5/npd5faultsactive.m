sector        = 'NPD5';
grdecl        = readGRDECL([sector, '.grdecl']);
actnum        = grdecl.ACTNUM;
grdecl.ACTNUM = ones(prod(grdecl.cartDims),1);
G             = processGRDECL(grdecl, 'checkgrid', false);

clf, subplot('position',[0.025 0.025 0.95 0.95]);
plotGrid(G,'FaceColor','none','EdgeAlpha',0.1);
plotFaces(G,find(G.faces.tag>0),'FaceColor','r');
axis tight off; view(-145,60);

plotGrid(G,find(actnum(G.cells.indexMap)), ...
         'FaceColor', 'b', 'FaceAlpha', 0.4, 'EdgeAlpha', 0.1);
view(20,75);