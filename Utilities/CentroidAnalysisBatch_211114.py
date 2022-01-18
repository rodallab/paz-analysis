# -*- coding: utf-8 -*-
"""
Created on Thu Aug 12 18:35:46 2021

@author: Steve
"""


# Iterate through csv files
# For each csv file, channel xy and NMJ xys
# create alpha mask with NMJxys
# For each channel and paired fake data, measure k and g
# aggregate across channels for each folder and create graph

import csv
import numpy as np
import os
from pathlib import Path
import matplotlib
import matplotlib.pyplot as plt
from pointpats import distance_statistics, PointPattern, PointProcess
import libpysal
import numba
import pointpats as pp
from descartes import PolygonPatch
import importlib
fn='Z:/Current members/DelSignore/PAZ Organization Analysis Project/PAZ/SD555_MaxIP-Control/test/AnalysisVSD555_V151/Centroids/SD555_5mm_L2__2020_12_06__10_33_28_Out-1_MaxZ14-Z24__Centroids.csv'
PATH='test'
COLORS=['g', 'm', 'k', 'c']
SPOTS=['.', '+', 'x']
SCALE=23.5641
DTRANSFORM=False   #whether to transform d values by multiplying by density prior to averaging
#CHANS=["Dyn", "BRP", "Nwk"]

class MaternPointProcess(PointProcess):
    """
    """

    
    def __init__(self, window, n, samples, keep=False,
                 asPP=False, conditioning=False):
        self.window = window
        self.n = n
        self.scale = 1 # scaling to create buffer to manage edge effects
        self.samples = samples
        self.parents = 2*n*scale*scale
        self.conditioning = conditioning
        
        lam = self.parents/window.area
        d = np.sqrt(1.594/lam/np.pi)
        self.radius = d
        self.keep = keep
        self.realizations = {}
        self.setup()
        for sample in range(samples):
            self.realizations[sample] = self.draw(self.parameters[sample])
        if asPP:
            for sample in self.realizations:
                points = self.realizations[sample]
                self.realizations[sample] = PointPattern(points, window=self.window)


    def setup(self):
        """
        Generate the number of events for each realization. If
        "conditioning" is False, all the event numbers are the same;
        if it is True, the number of parents is a random variable
        following a Poisson distribution, resulting in varied number
        of events.

        """

        self.parameters = {}
        self.num_parents = {}
        if self.conditioning:
            lambdas = poisson(self.parents, self.samples)
            for i, l in enumerate(lambdas):
                num = l * self.children
                self.parameters[i] = {'n': num}
                self.num_parents[i] = l
        else:
            for i in range(self.samples):
                self.parameters[i] = {'n': self.n}
                self.num_parents[i] = self.parents

    def draw(self, parameter):
        """
        Generate a series of point coordinates within the given window.

        Parameters
        ----------
        parameter  : dictionary
                     Key: 'n'.
                     Value: size of the realization.

        Returns
        -------
                   : array
                     A series of point coordinates.

        """
        c = 0
        sample = []
        n = parameter['n']
        while c < n:
            pnts = self.realize(n)
            pnts = [libpysal.cg.shapes.Point((x, y)) for x, y in pnts]
            pins = self.window.filter_contained(pnts)
            
            # repeat thinning on extended sample
            sample.extend(pins)
            revpnts = self.matII(np.array(sample).T)
            revpnts = [libpysal.cg.shapes.Point((x, y)) for x, y in revpnts]
            c = len(revpnts)
        return np.array([np.asarray(p) for p in revpnts[:n]])

    def realize(self, n):

        l, b, r, t = self.window.bbox
        x_buffer = int((r-l)*(self.scale-1)/2)
        y_buffer = int((t-b)*(self.scale-1)/2)
        nparents = int(2*n*scale*scale)
        
        # get parent points
        pxs = np.random.uniform(l-x_buffer, r+x_buffer, (nparents))
        pys = np.random.uniform(b-y_buffer, t+y_buffer, (nparents))
        pxy = np.array((pxs, pys))
        points = self.matII(pxy)
        return points
        
           
    def matII(self, points):
        # Thin according to Matern II process
        npoints = points[0].size
        ages = np.random.uniform(size=npoints)
        cents = np.array((points[0], points[1], ages))
        keep=np.zeros(npoints, dtype=bool)
        for ix in range(npoints):
            # return distance between each parent point ix and all other points
            distances=np.array([np.hypot(cents[0][ix]-cents[0][comp], cents[1][ix]-cents[1][comp]) for comp in range(npoints)])

            #check whether each point is inside the disc (but not the parent)
            inDisc=(distances<self.radius)&(distances>0)
            #get the age of the youngest point
            if(np.count_nonzero(inDisc)>0):
                youngest=np.min(cents[2][inDisc])
                if(cents[2][ix]<youngest):
                    keep[ix]=True
            else: keep[ix]=True
        keepers = cents[:2, keep]
        points=[(keepers[0][i], keepers[1][i]) for i in range(len(keepers[0]))]
        return points


def importData(fn):
    ''' 
    imports data from _centroids.csv file generated by PAZ analysis macro
    Expects 2*nChannels columns of local min/max centroid data followed
    by 2*nChannels columns of centroid data.
    Returns tuple of min/max and centroid coordinate lists
    '''
    mcent={}
    gcent={}
    with open(fn) as csvfile:
        csvr = csv.reader(csvfile, quotechar='|')
        header=next(csvr)
        # number of channels=(len(header))/4 
        # min/max and geometric centroid x and y for each channel  
        # last two columns are the NMJ coordinates
        table=np.array([row for row in csvr])
            
    
    chans=[header[i].split('_')[1] for i in range(0, int(len(header)/2), 2)]
    
    for c,chan in enumerate(chans):
        mcent[chan]=np.array([[row[2*c], row[2*c+1]] 
                              for row in table 
                              if row[2*c]!='NaN'], 
                             dtype=float)
        
        gcent[chan]=np.array([[row[2*len(chans)+2*c], row[2*len(chans)+2*c+1]] 
                              for row in table 
                              if row[2*len(chans)+2*c]!='NaN'], 
                             dtype=float)
    
    
    return (mcent, gcent)


def makeWindow(shape):
    alphaVerts=list(shape.exterior.coords)
    alphaVerts=[libpysal.cg.Point(p) for p in alphaVerts]
    alphaPoly=libpysal.cg.Polygon(alphaVerts)
    win=pp.Window(alphaPoly.parts)
    return win


def plotCentroids(centroids, mask=None, savepath=None):
    
        
    f,ax=plt.subplots(1,1)
    if(mask):
        ax.add_patch(
        PolygonPatch(
            mask,
            edgecolor='powderblue',
            facecolor='powderblue',
            alpha=.4,
            label='alpha shape'
        ))
    
    for i, (chan, cent) in enumerate(centroids.items()):
        col=i % len(COLORS)
        mark=i // len(COLORS)
        ax.scatter(cent[:,0],
                   cent[:,1], 
                   color=COLORS[col], 
                   marker=SPOTS[mark], 
                   label='{}'.format(chan)
                   )
        
    plt.legend()
   
    if(savepath):
        plt.savefig(os.path.join(savepath, "CentroidPositions.svg"), dpi=None)
    else:
        plt.show()
        
def plotG(obs, chan, pois=None, clust=None, spaced=None, savepath=None, flag=""):
    
    ax=plt.subplot()
    ax.plot(obs.d, obs.G, label = '{}'.format(chan), color='black', linewidth=2)
    
    if(pois):
        ax.plot(pois.d, pois.G, color='darkgray', linewidth=2, 
             label='Poisson')
    
    if(clust):
        ax.plot(clust.d, clust.G, color='magenta', linewidth=2, linestyle=(0, (5, 10)), label='Clustered')
        
    if(spaced):
        ax.plot(spaced.d, spaced.G, color='magenta', linewidth=2, label='Spaced')
    plt.legend()
    
    if(savepath):
        plt.savefig(os.path.join(savepath, "{}-{}_gd.svg".format(chan, flag)), dpi=None, facecolor='w', edgecolor='w')
    else:
        plt.show()


def runme(path):
    fl=os.listdir(path)
    allObsg={}
    allPoisg={}
    allClg={}
    allSpg={}
    allDens={}
    
    for fn in fl:
        print(fn)
        # Make an image specific folder to store some data
        imgfolder=os.path.join(path, Path(fn).stem)
        if(os.path.isdir(imgfolder)==False):
            os.mkdir(imgfolder)
            
        if(fn.endswith(".csv")):
            mcents, gcents = importData(os.path.join(path, fn))
            chans = mcents.keys()
            
            # Make sure each g function list has right number of slots
            if not allObsg:
                for c in chans:
                    allObsg[c]=[]
                    allPoisg[c]=[]
                    allClg[c]=[]
                    allSpg[c]=[]
                    allDens[c]=[]
                    
            # Make combined list to generate alpha mask & convert to window
            allcents=np.vstack([vals for vals in mcents.values()])
            alpha_shape=libpysal.cg.alpha_shape_auto(allcents)
            win=makeWindow(alpha_shape)
            
            # In img specific folder:
            # Make a plot of alpha mask and points
            # Make plot for each channel with model and observed g
            
            plotCentroids(mcents, alpha_shape, imgfolder)
            
            for i, (chan, centroids) in enumerate(mcents.items()):

                # Generate poisson, clustered, and spaced distributions
                area=alpha_shape.area
                dens=len(centroids)/area
                allDens[chan].append(dens)
                poisson=pp.PoissonPointProcess(win, len(centroids), 1, asPP=True)
                spaced = MaternPointProcess(win, len(centroids), 1, asPP=True)
                clustered = pp.PoissonClusterPointProcess(win, 
                                                          len(centroids), 
                                                          int(len(centroids)**.5), 
                                                          area/512, 1, 
                                                          asPP=True
                                                          )
                
                spcent=spaced.realizations[0]
                clcent=clustered.realizations[0]
                pscent=poisson.realizations[0]
                
                # Calculate g(d) for each distribution and plot
                obsg=distance_statistics.G(pp.PointPattern(centroids), 40)
                poisg=distance_statistics.G(pscent, 40)
                clg=distance_statistics.G(clcent, 40)
                spg=distance_statistics.G(spcent, 40)
                
                # Here multiply d by density. This transform seems to do a good job
                # of reducing variability between NMJs
                
                obsgd=obsg.d*dens
                poisd=poisg.d*dens
                cld=clg.d*dens
                spd=spg.d*dens
                plotG(obsg, chan, poisg, clg, spg, savepath=imgfolder)
                plt.clf()
                
                if(DTRANSFORM):
                    allObsg[chan].append((obsgd, obsg.G))
                    allPoisg[chan].append((poisd, poisg.G))
                    allClg[chan].append((cld, clg.G))
                    allSpg[chan].append((spd, spg.G))
                else:
                    allObsg[chan].append((obsg.d, obsg.G))
                    allPoisg[chan].append((poisg.d, poisg.G))
                    allClg[chan].append((clg.d, clg.G))
                    allSpg[chan].append((spg.d, spg.G))
                    
    meanObs, meanPois, meanClg, meanSpg, meanDens = {}, {}, {}, {}, {}
    
    for c in chans:
        # calculate the mean d and mean g for each set
        meanObs[c]=np.mean(allObsg[c], 0).transpose()
        meanPois[c]=np.mean(allPoisg[c], 0).transpose()
        meanClg[c]=np.mean(allClg[c], 0).transpose()
        meanSpg[c]=np.mean(allSpg[c], 0).transpose()
        meanDens[c]=np.mean(allDens[c])
        if DTRANSFORM:
            dcorr=meanDens[c]
        else:
            dcorr=1
        # store d and g as attributes to make compatible with plotG
        # restore units of d to pixels, and scale to microns
        meanObs[c] = type('obj', (object,), {'d' : meanObs[c][:,0]/dcorr/SCALE, 'G':meanObs[c][:,1]})
        meanPois[c] = type('obj', (object,), {'d' : meanPois[c][:,0]/dcorr/SCALE, 'G':meanPois[c][:,1]})
        meanClg[c] = type('obj', (object,), {'d' : meanClg[c][:,0]/dcorr/SCALE, 'G':meanClg[c][:,1]})
        meanSpg[c] = type('obj', (object,), {'d' : meanSpg[c][:,0]/dcorr/SCALE, 'G':meanSpg[c][:,1]})
        
        # Plot
        plotG(meanObs[c], 
              c, 
              meanPois[c], 
              meanClg[c], 
              meanSpg[c], 
              savepath=path, 
              flag="_{}-Avg".format(c)
              )
        plt.clf()
    

            