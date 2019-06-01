#############################################################################
#
# RandField2.py - a spatially-correlated random field generator in 2D or 3D
#
# by Walt McNab
#
#############################################################################


from numpy import *
from pandas import *
from scipy.interpolate import Rbf
from scipy.interpolate import NearestNDInterpolator
from scipy.spatial import distance
import scipy.stats as stats


class Params:
    
    def __init__(self):
        # miscellaneous setup parameters
        lineInput = []        
        inputFile = open('params.txt','r')
        for line in inputFile: lineInput.append(line.split())
        inputFile.close()
        self.grid0 = array([float(lineInput[1][1]), float(lineInput[1][2]), float(lineInput[1][3])])
        self.gridend = array([float(lineInput[2][1]), float(lineInput[2][2]), float(lineInput[2][3])])        
        self.n = array([int(lineInput[3][1]), int(lineInput[3][2]), int(lineInput[3][3])])
        self.a = array([float(lineInput[4][1]), float(lineInput[4][2]), float(lineInput[4][3])])
        self.depth = float(lineInput[5][1])
        self.stdev0 = float(lineInput[6][1])
        self.lower = float(lineInput[7][1])
        self.upper = float(lineInput[8][1])
        self.rsearch0 = float(lineInput[9][1])
        self.expn = float(lineInput[10][1])        
        self.dmin = float(lineInput[11][1])
        self.f = lineInput[12][1]
        print('Read setup parameters.')


class Grid:
    
    def __init__(self, params):

        # seed points
        points = read_csv('seeds.csv', sep=',')
        x = array(points['x'])
        y = array(points['y'])
        z = array(points['z'])
        v = array(points['v'])        
        
        # grid setup
        self.lengthScale = array([params.gridend[0]-params.grid0[0],
            params.gridend[1]-params.grid0[1],
            params.gridend[2]-params.grid0[2]])
        self.dcell = self.lengthScale/params.n
        xgrid = linspace(params.grid0[0]+0.5*self.dcell[0], params.gridend[0]-0.5*self.dcell[0], params.n[0])
        ygrid = linspace(params.grid0[1]+0.5*self.dcell[1], params.gridend[1]-0.5*self.dcell[1], params.n[1])
        zgrid = linspace(params.grid0[2]+0.5*self.dcell[2], params.gridend[2]-0.5*self.dcell[2], params.n[2])        
        Z, Y, X = meshgrid(zgrid, ygrid, xgrid, indexing='ij')
        xg = X.flatten()
        yg = Y.flatten()
        zg = Z.flatten()
        marked = zeros(len(xg), bool) * False
        val = zeros(len(xg), float) - 9999.
        self.cells = DataFrame(data={'x':xg, 'y':yg, 'z':zg, 'v':val, 'marked':marked})
        print('Read seed points and set up grid data frame.')

        # populate grid cells that contain the initial seed points
        print('Marking grid cells containing seed points ...')
        indx = self.GetIndex(x, y, z, params)
        self.cells['marked'].iloc[indx] = True
        self.cells['v'].iloc[indx] = v
        
    def GetIndex(self, x, y, z, params):
        # index number of grid cell corresponding to (x, y, z)
        col = ((x-params.grid0[0])/self.dcell[0]).astype(int)
        row = ((y-params.grid0[1])/self.dcell[1]).astype(int)
        layer = ((z-params.grid0[2])/self.dcell[2]).astype(int)
        indx = col + row*params.n[0] + layer*params.n[0]*params.n[1]
        return indx

    def Interp(self, x, y, z, v, xp, yp, zp, params, d):
        # search range-dependent estimate for v at (xp, yp, zp)
        rsearch = params.rsearch0 * d     # reduce size of subset to reflect updated (effective) search radius
        tpts = transpose([x, y, z])
        dist = distance.cdist([[xp, yp, zp]], tpts)
        x = x[dist[0]<=rsearch]
        y = y[dist[0]<=rsearch]
        z = z[dist[0]<=rsearch]
        v = v[dist[0]<=rsearch]
        rbfi = Rbf(x, y, z, v, function=params.f)   # radial basis function interpolator
        return rbfi(xp, yp, zp)

    def Fill(self, params):
        # fill in un-marked grid cells by interpolation
        print('Filling remaining cells by interpolation ...')
        marked = self.cells[self.cells['marked']==True].copy()
        marked.to_csv('points.csv', index=False)
        x = array(marked['x']) / params.a[0]
        y = array(marked['y']) / params.a[1]        
        z = array(marked['z']) / params.a[2]
        v = array(marked['v'])
        fN = NearestNDInterpolator((x, y, z), v)
        unmarked = self.cells[self.cells['marked']==False].copy()
        xp = array(unmarked['x']) / params.a[0]
        yp = array(unmarked['y']) / params.a[1]
        zp = array(unmarked['z']) / params.a[2]
        vp = fN(xp, yp, zp)
        unmarked['marked'] = True
        unmarked['v'] = vp
        self.cells = concat([marked, unmarked])

    def Posit(self, params, d):
        # return a posited new data point (represented by grid cell)
        subset = self.cells[self.cells['marked']==False]
        indices = list(subset.index.values)
        r = random.choice(indices)          # select a random grid cell to assign value
        rcell = subset.loc[r]
        xp = rcell['x'] / params.a[0]
        yp = rcell['y'] / params.a[1]
        zp = rcell['z'] / params.a[2]
        marked = self.cells[self.cells['marked']==True]
        x = array(marked['x']) / params.a[0]
        y = array(marked['y']) / params.a[1] 
        z = array(marked['z']) / params.a[2]
        v0 = array(marked['v'])
        mu = self.Interp(x, y, z, v0, xp, yp, zp, params, d)
        sigma = params.stdev0 * d
        X = stats.truncnorm((params.lower - mu) / sigma, (params.upper - mu) / sigma, loc=mu, scale=sigma)
        v = X.rvs(1)[0]
        self.cells['marked'].iloc[r] = True
        self.cells['v'].iloc[r] = v        

    def WriteOutput(self):
        # output to file
        self.cells = self.cells[['x', 'y', 'z', 'v']]
        self.cells.to_csv('filled_cells.csv', index=False)


def RandField():

    # read model parameters
    params = Params()
    
    # read seed points and construct grid
    grid = Grid(params)    
 
    # step through cells
    print('Positing pilot points ...')
    nfilled = int(params.depth*prod(params.n))
    for i in range(nfilled):
        d = 1.0 - (1.0-params.dmin)*(float(i)/float(nfilled))**params.expn
        print('count = ' + str(i) + '/' + str(nfilled), 'search = ' + str(d))          
        grid.Posit(params, d)
       
    # fill in remaining cells with straight interpolation of existing marked cell set
    grid.Fill(params)
    
    # write completed point set
    grid.WriteOutput()
    
    print('Done.')


### run script ###
RandField()
