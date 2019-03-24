sector        = 'NPD5';
grdecl        = readGRDECL([sector, '.grdecl']);
actnum        = grdecl.ACTNUM;
grdecl.ACTNUM = ones(prod(grdecl.cartDims),1);
G             = processGRDECL(grdecl, 'checkgrid', false);

clf
K = reshape(load([sector, '_Permeability.txt']')', prod(G.cartDims), []);
rock.perm = bsxfun(@times, [1 1 0.1], K(G.cells.indexMap)).*milli*darcy; clear K;
hp = plotCellData(G,log10(rock.perm(:,1)),'EdgeColor','k','EdgeAlpha',0.1);
view(-45,15), axis tight off, zoom(1.2)

% Manipulate the colorbar to get the ticks we want
h = colorbar;
cs = [0.01 0.1 1 10 100 1000];
caxis(log10([min(cs) max(cs)]*milli*darcy));
set(h, 'XTick', 0.5, 'XTickLabel','mD', ...
   'YTick', log10(cs*milli*darcy), 'YTickLabel', num2str(cs'));

delete(hp), view(-20,35)
plotGrid(G,'FaceColor','none','EdgeAlpha',0.1);
hp = plotCellData(G,log10(rock.perm(:,1)), ...
                  find(rock.perm(:,1)>0.01*milli*darcy), ...
                  'EdgeColor','k', 'EdgeAlpha', 0.1);