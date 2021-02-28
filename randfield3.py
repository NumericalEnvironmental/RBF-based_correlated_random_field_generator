#############################################################################
#
# RandField3.py - a spatially-correlated random field generator in 2D or 3D
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

class InvDist:
    
    def __init__(self, xp, yp, zp, vp):
        self.pts = transpose([xp, yp, zp])
        self.vp = vp
        
    def Interpolate(self, location, params):   
        d0 = 1e-6
        d = distance.cdist(location, self.pts)[0] + d0
        pointSet = DataFrame(data={'distance':d, 'value':self.vp})
        nearSet = pointSet[pointSet['distance']<=params.searchInvDst]
        dLocal = sqrt(array(nearSet['distance']**2) + params.smooth**2)
        vLocal = array(nearSet['value'])
        h = sum(vLocal/dLocal**2) / sum(1./dLocal**2)
        return h


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
        self.mu = float(lineInput[6][1])
        self.stdev0 = float(lineInput[7][1])
        self.lower = float(lineInput[8][1])
        self.upper = float(lineInput[9][1])
        self.rsearch0 = float(lineInput[10][1])
        self.expn = float(lineInput[11][1])        
        self.dmin = float(lineInput[12][1])
        self.f = lineInput[13][1]
        self.searchInvDst = float(lineInput[14][1])
        self.smooth = float(lineInput[15][1])
        self.rescale = bool(lineInput[16][1])
        self.rescaleMu = float(lineInput[17][1])
        self.rescaleSigma = float(lineInput[18][1])        
        print('Read setup parameters.')


class CoVary:
    
    def __init__(self):
        # define variable that co-varies with main spatially-correlated variable
        lineInput = []        
        inputFile = open('covary.txt','r')
        for line in inputFile: lineInput.append(line.split())
        inputFile.close()
        self.b = float(lineInput[0][1])
        self.m = float(lineInput[1][1])        
        self.sigma = float(lineInput[2][1])
        self.fraction = float(lineInput[3][1])        
        print('Read co-variable parameters.')


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
        covar0 = empty(len(xg))
        covar0[:] = nan
        self.cells['cov'] = covar0      # placeholder for co-variable
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
        
    def Fill(self, mode, params):
        # fill in un-marked grid cells by interpolation
        if mode == 0:
            print('Filling remaining primary variable cells by interpolation ...')
            marked = self.cells[self.cells['marked']==True].copy()
            marked.to_csv('points.csv', index=False)    # write points to file
            v = array(marked['v'])
        else:
            print('Filling remaining co-variable cells by interpolation ...')
            marked = self.cells[~self.cells['cov'].isnull()].copy()
            v = array(marked['cov'])
        x = array(marked['x']) / params.a[0]
        y = array(marked['y']) / params.a[1]        
        z = array(marked['z']) / params.a[2]
        fN = InvDist(x, y, z, v)
        xp = array(self.cells['x']) / params.a[0]
        yp = array(self.cells['y']) / params.a[1]
        zp = array(self.cells['z']) / params.a[2]
        vp = []
        for i in range(len(xp)): vp.append(fN.Interpolate( [[xp[i], yp[i], zp[i]]], params))
        if mode == 0: self.cells['v'] = vp
        else: self.cells['cov'] = vp        

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
        x, y, z, v = DistFilter(x, y, z, v0, xp, yp, zp, params, d)  # find nearby points
        if len(x) > 1:
            rbfi = Rbf(x, y, z, v, function=params.f)   # radial basis function interpolator
            mu = rbfi(xp, yp, zp)
        else:
            mu = params.mu
        sigma = params.stdev0 * d
        X = stats.truncnorm((params.lower - mu) / sigma, (params.upper - mu) / sigma, loc=mu, scale=sigma)
        v = X.rvs(1)[0]
        self.cells['marked'].iloc[r] = True
        self.cells['v'].iloc[r] = v        

    def CovNoise(self, covary):
        # populate co-variable in selected subset of marked points
        marked = self.cells[self.cells['marked']==True].copy()
        indices = list(marked.index.values)
        nSet = int(covary.fraction * len(marked))
        rSet = random.choice(indices, size=nSet)
        rCell = marked.loc[rSet]
        noise = random.normal(0., covary.sigma, len(rSet))
        cv = covary.b + covary.m*rCell['v'] + noise
        R2 = corrcoef(array(rCell['v']), cv)[0, 1]**2
        print('Point correlation (R^2) =', R2)
        self.cells['cov'].iloc[rSet] = cv

    def WriteOutput(self):
        # output to file
        self.cells = self.cells[['x', 'y', 'z', 'v', 'cov']]
        self.cells.to_csv('filled_cells.csv', index=False)


### utility functions ###


def DistFilter(x, y, z, v, xp, yp, zp, params, d):
    # filter data set by distance
    rsearch = params.rsearch0 * d     # reduce size of subset to reflect updated (effective) search radius
    tpts = transpose([x, y, z])
    dist = distance.cdist([[xp, yp, zp]], tpts)
    x = x[dist[0]<=rsearch]
    y = y[dist[0]<=rsearch]
    z = z[dist[0]<=rsearch]
    v = v[dist[0]<=rsearch]
    return x, y, z, v


### main script ###


def RandField():

    # read model parameters
    params = Params()
    
    # read seed points and construct grid
    grid = Grid(params)    
 
    # step through cells
    print('Positing pilot points ...')
    nfilled = int(params.depth*params.n[0]*params.n[1]*params.n[2])
    for i in range(nfilled):
        d = 1.0 - (1.0-params.dmin)*(float(i)/float(nfilled))**params.expn
        print('count = ' + str(i) + '/' + str(nfilled), 'search = ' + str(d))          
        grid.Posit(params, d)
       
    # compute points for co-variable
    covary = CoVary()
    grid.CovNoise(covary)
       
    # fill in remaining cells with straight interpolations of existing marked cell set
    grid.Fill(0, params)   # primary variable
    
    # re-scale primary variable by stretching histogram (assumes a normal distribution)
    if params.rescale:
        u0 = grid.cells['v'].mean()
        stdev0 = grid.cells['v'].std()
        vcdf = stats.norm.cdf(grid.cells['v'], loc=u0, scale=stdev0)
        vScaled = stats.norm.ppf(vcdf, loc=params.rescaleMu, scale=params.rescaleSigma)
        shift = vScaled - grid.cells['v']
        grid.cells['v'] = vScaled
    
    # process covariable
    grid.Fill(1, params)   # co-variable
    if params.rescale: grid.cells['cov'] = grid.cells['cov'] + shift*covary.m       # re-scale co-variable
    
    # write completed point set
    grid.WriteOutput()
    
    print('Done.')


### run script ###
RandField()
