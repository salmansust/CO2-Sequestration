sector        = 'NPD5';
grdecl        = readGRDECL([sector, '.grdecl']);
actnum        = grdecl.ACTNUM;
grdecl.ACTNUM = ones(prod(grdecl.cartDims),1);
G             = processGRDECL(grdecl, 'checkgrid', false);

w = load([sector, '_Well.txt']);
W = verticalWell([], G, rock,  w(1,1), w(1,2), w(1,3):w(1,4),  ...
                 'Radius', 0.1, 'name', 'I');
plotWell(G,W,'height',1000,'color','r');
zoom out