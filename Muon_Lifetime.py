import numpy as np
import matplotlib.pyplot as plt
import pylab as pl
import random
from tqdm import tqdm
import math
from scipy import stats
from __future__ import print_function

### Function that finds intersection between line segment and plane
def LinePlaneCollision(planeNormal, planePoint, rayDirection, rayPoint, epsilon=1e-6):
 
	ndotu = planeNormal.dot(rayDirection)
	if abs(ndotu) < epsilon:
		raise RuntimeError("no intersection or line is within plane")
 
	w = rayPoint - planePoint
	si = -planeNormal.dot(w) / ndotu
	Psi = w + si * rayDirection + planePoint
	return Psi

# Varies size of the muon emission plane
heightplanelist = np.arange(1200,0,-100)
heightfluxvalues = []

for heightplane in tqdm(heightplanelist):
    xmin = -np.pi/2
    xmax = np.pi/2
    N = 2000000 # 2 million
  
    # Defines probability distribution function (pdf) that the Von Neumann Monte Carlo generation will follow
    def pdf(range):
        pdf_array = (np.cos(range))**2
        return pdf_array

    x=np.linspace(xmin,xmax,1000)

    pmin=0
    pmax=1

    naccept=0  
    ntrial=0

    # Von Neumann Monte Carlo generation
    ran=[]
    while naccept<N:
    
        x=np.random.uniform(xmin,xmax) # x'  
        y=np.random.uniform(pmin,pmax) # y'
    
        if y<pdf(x):
            ran.append(x)  
            naccept=naccept+1
        ntrial=ntrial+1
        
    ran=np.asarray(ran)
    
    ### Constraints for plane where muons are emitted
    ###
    size = 2000 # 20 m
    xplanemin = -size/2
    xplanemax = size/2
    yplanemin = -size/2
    yplanemax = size/2
    height = 12.5 # height of detector
    radius = 7.5 # radius of detector
    ###

    ### Generating N points on plane for muon to be emitted at
    xyplane = []
    for i in range(0,N):
        xyplane.append([np.random.uniform(xplanemin,xplanemax),np.random.uniform(yplanemin,yplanemax),heightplane])
        
    ### Muon momentum, need to calculate unit vector of muon emission
    ### r = 1
    ### x = sin(theta) cos(phi)
    ### y = sin(theta) sin(phi)
    ### z = cos(theta)
    ptheta = ran

    # Generating random phi momenta for each theta momentum value
    pphi = []
    for i in range(0,N):
        pphi.append(np.random.uniform(0,2*np.pi))

    # Now the unit vector can be calculated in cartesian coordinates
    px = np.sin(ptheta) * np.cos(pphi)
    py = np.sin(ptheta) * np.sin(pphi)
    pz = -np.cos(ptheta)

    # Creates nested array containing unit vectors
    unitv = []
    for i in range(0,N):
        unitv.append([px[i],py[i],pz[i]])
        
    ### Number of planes we are considering from z=0 -> z=h (h is height of detector)
    ### Code that selects optimum number of planes so that all muons will be detected
    ###
    for n in range(2,1001):
        thetar = math.atan((heightplane - (n-1)/n * height)/(xplanemax + radius))

        thetad = thetar*180/np.pi

        phir = math.atan(height/(2*radius*n))

        phid = phir*180/np.pi

        if phid<thetad:
            planenumber = n
            break
    ###
    planeheight = np.arange(height,0-(height/planenumber),-(height/planenumber))

    # Sets initial flux values to zero
    topflux = 0
    sideflux = 0

    planepointarray = []

    for i in planeheight: # For each plane
        # Define plane
        planenormal = np.array([0, 0, 1])
        planepoint = np.array([0, 0, i]) #Any point on the plane
        planepointarray.append(planepoint)
        unitv_temp=[]
        xyplane_temp=[]

        for j in range(0,len(unitv)): # For each ray

            # Define ray
            raydirection = np.array(unitv[j])
            raypoint = np.array(xyplane[j])

            # Calculates intercept between muon and xy plane at selected height
            intercept = LinePlaneCollision(planenormal, planepoint, raydirection, raypoint)

            # Determines whether the muon has hit the inside of the detector (or top for first plane)
            if intercept[0]**2 + intercept[1]**2 <= radius**2:

                # Appends flux values
                if i == 12.5:
                    topflux = topflux + 1
                else:
                    sideflux = sideflux + 1

            # Adds unit vectors and plane point pairs for muons that didn't hit the detector at the
            # selected height and removes muon that did hit it to not double count.
            else:
                unitv_temp.append(unitv[j])
                xyplane_temp.append(xyplane[j])

        # Updates next array of unit vector and plane point pairs for the next xy plane
        unitv = unitv_temp
        xyplane = xyplane_temp
        
    heightfluxvalues.append([topflux,sideflux])
    print("For plane height =",heightplane/100,"m , number of muons =",N)
    print([topflux,sideflux])
