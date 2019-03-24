sector        = 'NPD5';
grdecl        = readGRDECL([sector, '.grdecl']);
actnum        = grdecl.ACTNUM;
grdecl.ACTNUM = ones(prod(grdecl.cartDims),1);
G             = processGRDECL(grdecl, 'checkgrid', false);

clf
p = reshape(load([sector, '_Porosity.txt'])', prod(G.cartDims), []);
rock.poro = p(G.cells.indexMap); clear p
hp = plotCellData(G,rock.poro,'EdgeColor','k','EdgeAlpha',0.1);
colorbar; caxis([0.1 0.3]), view(-45,15), axis tight off, zoom(1.2)

delete(hp), view(-15,40)
plotGrid(G,'FaceColor','none','EdgeAlpha',0.1);
plotCellData(G,rock.poro, find(rock.poro>0.1), ...
             'EdgeColor','k','EdgeAlpha',0.1);