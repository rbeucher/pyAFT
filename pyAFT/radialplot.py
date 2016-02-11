import pylab as plt
import random
import numpy as np
from matplotlib import collections  as mc
from matplotlib.ticker import MaxNLocator

def _radius(a, se, f):
    """Return Ellipse radius at max precision position"""
    r = max(np.sqrt(se**2 + (f*a)**2))
    return r

def _f(xlim, ylim):
    # Calculate conversion factor for plot coordinates
    return (max(xlim) - min(xlim)) / (max(ylim)- min(ylim))

def _zlim(values):
    zspan = (np.mean(values) *0.5) / (np.std(values) * 100)
    zspan = 0.9 if zspan > 1.0 else zspan
    zlim = (np.floor((0.9 - zspan)*min(values)),
            np.ceil((1.1 + zspan)*max(values)))
    return zlim

def plot_line(val, central, r, f, **kwargs):
    ax = plt.gca()
    x1,y1 = (0,0)
    x2 = r / np.sqrt(1 + f**2*(val-central)**2)
    y2 = x2*(val-central)
    ax.plot((x1,x2), (y1, y2), **kwargs)
    
def get_ellipse_coord(val,r,f,central):
    x = r / np.sqrt(1 + f**2*(val - central)**2)
    y = x * (val - central)
    return x, y

def plot_ticks(x1,y1,x2,y2):
    ax = plt.gca()
    A = list(zip(x1, y1))
    B = list(zip(x2, y2))
    segments = list(zip(A, B))
    lc = mc.LineCollection(segments)
    ax.add_collection(lc)
    lc.set_color("k")



class Radialplot():

    def __init__(self, values, errors, central):

        self.values = values
        self.errors = errors
        self.central = central

        # Calculate standard estimates and z-values
        self.se = 1.0 / errors
        self.z = (values - central) / errors
        
        # Define Axes limits
        self.zticks_major = 1.015
        self.zticks_minor = 1.007
        self.zlabels = 1.02
        self.rellipse = 1.1
        
        minx = 0.0
        maxx = 1.2*max(self.se)
        miny = -15
        maxy = 15
        
        self.xlim=(minx,maxx)
        self.ylim=(miny,maxy)
        
        self.f = _f(self.xlim, self.ylim)
        # Now we need to create a z-axis
        self.r = self.rellipse*_radius(self.z, self.se, self.f)
        # Calculate z-span
        self.zlim = _zlim(values)
        
        
        # Create z-ticks
        locator = MaxNLocator(integer="true", nbins=5, prune="both", symetric="True")
        self.ticks_values_major = locator.tick_values(*self.zlim)
        locator = MaxNLocator(integer="true", nbins=20, prune="both", symetric="True")
        self.ticks_values_minor = locator.tick_values(*self.zlim)
        
        # Calculate major z-ticks coordinates
        self.ticks_x1_major, self.ticks_y1_major = get_ellipse_coord(self.ticks_values_major, 
                                                                     self.r, self.f, 
                                                                     self.central)
        self.ticks_x2_major, self.ticks_y2_major = get_ellipse_coord(self.ticks_values_major, 
                                                                     self.zticks_major*self.r,
                                                                     self.f,
                                                                     self.central)
        
        # Calculate minor z-ticks coordinates
        self.ticks_x1_minor, self.ticks_y1_minor = get_ellipse_coord(self.ticks_values_minor, 
                                                                     self.r,self.f,
                                                                     self.central)
        self.ticks_x2_minor, self.ticks_y2_minor = get_ellipse_coord(self.ticks_values_minor, 
                                                                     self.zticks_minor*self.r,
                                                                     self.f,
                                                                     self.central)
        
        # Calculate z-labels positions
        self.xlab, self.ylab = get_ellipse_coord(self.ticks_values_major,self.zlabels*self.r,
                                                 self.f, self.central)
        
        # Ellipse
        ellipse_values = np.linspace(min(np.hstack((self.ticks_values_major,
                                                    self.ticks_values_minor))), 
                                    max(np.hstack((self.ticks_values_major, 
                                                   self.ticks_values_minor))), 500)
        
        self.ellipse_x, self.ellipse_y =  get_ellipse_coord(ellipse_values,
                                                            self.r, self.f, 
                                                            self.central)
        
   
    def plot(self):
        ## Plotting
        fig, ax = plt.subplots()
        
        
        ax.set_ylim(self.ylim)
        ax.set_xlim(self.xlim)
        ax.set_yticks([-2,0,2])
        ax.spines["left"].set_bounds(-2,2)
        ax.spines["bottom"].set_bounds(self.xlim[0], max(self.se))
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        # Only show ticks on the left and bottom spines
        ax.yaxis.set_ticks_position('left')
        ax.xaxis.set_ticks_position('bottom')
        ax.xaxis.set_tick_params(direction="in", pad=-15)
        
        # Now create a second X axis
        #newax = ax.twiny()
        #newax.set_xlim(xlim)
        #newax.set_ylim(ylim)
        #newax.set_frame_on(True)
        #newax.patch.set_visible(False)
        #newax.xaxis.set_ticks_position('bottom')
        #newax.xaxis.set_label_position('bottom')
        #newax.spines["bottom"].set_bounds(xlim[0], max(se))
        #newax.set_aspect(f)
        # Only show ticks on the left and bottom spines
        #newax.xaxis.set_ticks_position('bottom')
        #newax.xaxis.set_tick_params(direction="out")
                
        ax.plot(self.ellipse_x, self.ellipse_y, c="k")
        plot_ticks(self.ticks_x1_major, self.ticks_y1_major,
                   self.ticks_x2_major, self.ticks_y2_major)
        plot_ticks(self.ticks_x1_minor, self.ticks_y1_minor,
                   self.ticks_x2_minor, self.ticks_y2_minor)
        
        
        # Add labels z-axis
        for index, label in enumerate(zip(self.xlab, self.ylab)):
            ax.text(*label, "{}".format(self.ticks_values_major[index]))
        
        # Plot data
        ax.plot(self.se, self.z, marker="o", linestyle="")
        
        # Plot central value line
        plot_line(self.central, self.central, self.r, self.f, c="k")

        ax.set_aspect(self.f)


if __name__== "__main__":

    # Create some dummy ages that follows a poisson distribution
    central = 200.
    values = np.random.normal(central,10, size=1000)
    errors = np.random.normal(10,1, size=1000)

    radialplot = Radialplot(values, errors, central)
    radialplot.plot()
    plt.show()
