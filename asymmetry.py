import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as col
import astropy.io.fits as pf
from astropy.io import ascii
import glob
import os
import operator
import _pickle as pickle
from scipy.ndimage import filters
from scipy.interpolate import interp1d
from scipy import stats
from matplotlib.backends.backend_pdf import PdfPages
import sep



def open_img(filename): 
	hdulist=pf.open(filename) 
	img=hdulist[0].data  
	return img

def exp_time(filename):
	hdulist=pf.open(filename) 
	hdr=hdulist[0].header 
	return hdr['EXPTIME']

def get_name(filename):
	i = 0
	while filename[i] != '_':
		i+=1
	obj_name = filename[20:i]
	return obj_name

def circle(radius,x0,y0):
	x = np.arange(100)
	y = x.reshape(-1,1)
	d = np.sqrt((x-x0)**2+(y-y0)**2)
	mask = d<radius
	return mask


data = '/home/conor/Galfit/fit_results_f140.txt'
data1 = '/home/conor/Galfit/fit_results_2f140.txt'
centroids = '/home/conor/Catalogs/q2343_centroids.txt'
datadir = '/home/conor/Project/'
ogimg = '/home/conor/Downloads/q2343_f140w_sci.fits'

fits_list = sorted(glob.glob(datadir+'*.fits'))





lines = ascii.read(data)
lines1 = ascii.read(data1)
lines2 = ascii.read(centroids)

gcx = np.array([])
gcy = np.array([])

x = np.array([])
y = np.array([])
cx = np.array([])
cy = np.array([])
objects = np.array([])
offset = np.array([])    # Distance from centroiding centroids from maxima
goffset = np.array([])   # Distance from galfit centroids from centroiding centroids
r70 = np.array([])
dc = np.array([])
ba = np.array([])

names = np.array([])
bovera = np.array([])

t1 = ['NB3024','NB2957','NB2910','NB2675','NB1684','NB1494','NB1041','NB0577','NB0193'] #I categorized the galaxies based on type, t2 has two galaxies close togehter while t1 appears to be disturbed
t2 = ['NB2929','NB2807','NB2571','NB2174','NB1579','NB1421','NB1420','NB1416','NB0970','NB0787','NB0743','NB0220']

for i in range(len(lines['Object'])):
	objects = np.append(objects,lines['Object'][i])
	gcx = np.append(gcx,lines['Centroid_X'][i])
	gcy = np.append(gcy,lines['Centroid_Y'][i])
	ba = np.append(ba,1 - lines['b/a'][i])



for i in range(len(lines1['Object'])):
	objects = np.append(objects,lines1['Object'][i])
	gcx = np.append(gcx,lines1['Centroid_X'][i])
	gcy = np.append(gcy,lines1['Centroid_Y'][i])
	ba = np.append(ba,1 - lines['b/a'][i])
	

for i in range(len(lines2['Object'])):
	obj = lines2['Object'][i]
	if obj not in objects:
		continue
	else:
		offset = np.append(offset,lines2['Offset'][i])
		x = np.append(x,lines2['X'])
		y = np.append(y,lines2['Y'])
		cx = np.append(cx,lines2['Centroid_X'][i])
		cy = np.append(cy,lines2['Centroid_Y'][i])
		r70 = np.append(r70,lines2['r70 est'][i])
		dc = np.append(dc,lines2['Distance'])


cxoffset = gcx - cx

cyoffset = gcy - cy

for i in range(len(cxoffset)):
	goffset = np.append(goffset,np.sqrt((cxoffset[i]**2 + cyoffset[i]**2)))


asymmetries = []

sizes = []



gasym = np.array([])

backgrounds = []
test1 = []

img = open_img(ogimg)
expt = exp_time(ogimg)





for filename in fits_list:
	img = open_img(filename)
	expt = exp_time(filename)

	rotated = np.rot90(img,2)

	

	name = get_name(filename)


	if name not in objects:
		continue


	i = np.where(objects == name)
	
	
	obj1 = objects[i]
	centroidx1 = cx[i]
	centroidy1 = cy[i]
	x1 = x[i]
	y1 = y[i]
	r70_1 = r70[i]

	print(name,centroidx1,centroidy1)

	asym1 = []

	

	names = np.append(names,name)
	bovera = np.append(bovera,ba[i])

	if name not in t1:
		r = 2*int(r70_1)

		raw = img[int(centroidx1)-r:int(centroidx1)+r,int(centroidy1)-r:int(centroidy1)+r]

		sizes.append(np.shape(raw)[0]**2)

		xmax = int(np.where(raw==np.amax(raw))[0][0]) + int(centroidx1) - r
		ymax = int(np.where(raw==np.amax(raw))[1][0]) + int(centroidy1) - r


		test = []


		positions = np.array([[],[]])

		for i in range(xmax+r+3,xmax+r+7):
			for j in range(ymax+r+3,ymax +r+7):
			

				raw = img[i-r:i+r,j-r:j+r]
				rotate = np.rot90(raw,1)
		        

				sub = np.abs(raw - rotate)
				resid = np.sum(sub)

					

				test.append(resid)



		background = np.min(test)



		for i in range(-2,2+1,1):
			for j in range(-2,2+1,1):
				

				

				raw = img[xmax+i-r:xmax+i+r,ymax+j-r:ymax+j+r]

				rotate = np.rot90(raw,1)

				flux = np.sum(np.abs(raw))
				if flux == 0:
					continue

				sub = np.abs(raw - rotate)

				resid = np.sum(sub)


				asym1.append(resid/flux-background/flux)

				


		if len(asym1) == 0:
			offset = np.delete(offset,i)
			names = np.delete(names,-1)
			continue



		asym = np.min(np.abs(asym1))
		
		
		r = 2
		index = np.where(np.abs(asym1)==asym)[0]


		index1 = int(index % (2*r+1))

		index2 = int((index - index1)/(2*r+1))


		i = index2 - r

		j = index1 - r

		r = 2*int(r70_1)

		raw = img[xmax+i-r:xmax+i+r,ymax+j-r:ymax+j+r]
		
		rotate = np.rot90(raw,1)

					
		flux = np.sum(np.abs(raw))
		if flux == 0:
			continue

		sub = np.abs(raw - rotate)
		resid = np.sum(sub)


		fig, (ax1,ax2,ax3) = plt.subplots(nrows=1,ncols=3)
		min = np.amin(raw)
		max = np.amax(raw)
		ax1.imshow(raw,cmap = 'binary',origin='lower',vmin = min, vmax = max)
		ax2.imshow(rotate,cmap='binary',origin='lower',vmin = min, vmax = max)
		m1 = ax3.imshow(sub,cmap='binary',origin='lower',vmin = min, vmax = max)
		cbar = fig.colorbar(m1, ax=ax3,shrink = .3)
		plt.text(0, -10, 'A = ' + str(asym))

		pp = PdfPages('/home/conor/Lyman/Asymmetries/model_'+name+'.pdf')
		pp.savefig(fig)
		pp.close()

		plt.close()
			

		gasym=np.append(gasym,asym)


		asymmetries.append(asym)
	if name in t1:
		r = int(10)

		raw = img[int(centroidx1)-r:int(centroidx1)+r,int(centroidy1)-r:int(centroidy1)+r]

		sizes.append(np.shape(raw)[0]**2)

		xmax = int(np.where(raw==np.amax(raw))[0][0]) + int(centroidx1) - r
		ymax = int(np.where(raw==np.amax(raw))[1][0]) + int(centroidy1) - r

		test = []


		positions = np.array([[],[]])

		for i in range(xmax+r+3,xmax+r+7):
			for j in range(ymax+r+3,ymax +r+7):
			

				raw = img[i-r:i+r,j-r:j+r]
				rotate = np.rot90(raw,1)
		        

				sub = np.abs(raw - rotate)
				resid = np.sum(sub)
					

				test.append(resid)



		background1 = np.min(test)



		for i in range(-3,3+1,1):
			for j in range(-3,3+1,1):
				

				raw = img[xmax+i-r:xmax+i+r,ymax+j-r:ymax+j+r]
				#rotate = rotated[int(centroidx1)+i-2*int(r70_1):int(centroidx1)+i+2*int(r70_1),int(centroidy1)+i-2*int(r70_1):int(centroidy1)+i+2*int(r70_1)]
				rotate = np.rot90(raw,1)

	
			
				flux = np.sum(np.abs(raw))
				if flux == 0:
					continue

				sub = np.abs(raw - rotate)


				'''fig, (ax1,ax2,ax3) = plt.subplots(nrows=1,ncols=3)
				ax1.imshow(raw,origin='lower')
				ax2.imshow(rotate,origin='lower')
				ax3.imshow(sub,origin='lower')
				plt.show()'''

				resid = np.sum(sub)

				

				asym1.append(resid/flux-background1/flux)


		if len(asym1) == 0:
			offset = np.delete(offset,i)
			names = np.delete(names,-1)
			continue



		asym = np.min(np.abs(asym1))

		

		
		
		r = 3
		index = np.where(np.abs(asym1)==asym)[0]


		index1 = int(index % (2*r+1))

		index2 = int((index - index1)/(2*r+1))


		i = index2 - r

		j = index1 - r

		r = 10

		raw = img[xmax+i-r:xmax+i+r,ymax+j-r:ymax+j+r]
		
		rotate = np.rot90(raw,1)

					
		flux = np.sum(np.abs(raw))
		if flux == 0:
			continue

		sub = np.abs(raw - rotate)
		resid = np.sum(sub)


		fig, (ax1,ax2,ax3) = plt.subplots(nrows=1,ncols=3)
		min = np.amin(raw)
		max = np.amax(raw)
		ax1.imshow(raw,cmap = 'binary',origin='lower',vmin = min, vmax = max)
		ax2.imshow(rotate,cmap='binary',origin='lower',vmin = min, vmax = max)
		m1 = ax3.imshow(sub,cmap='binary',origin='lower',vmin = min, vmax = max)
		cbar = fig.colorbar(m1, ax=ax3,shrink = .3)
		plt.text(0, -10, 'A = ' + str(asym))

		pp = PdfPages('/home/conor/Lyman/Asymmetries/model_'+name+'.pdf')
		pp.savefig(fig)
		pp.close()

		plt.close()
			

		gasym=np.append(gasym,asym)


		
		asymmetries.append(asym)



data = [names,asymmetries]


ascii.write(data,'/home/conor/Galfit/asymmetries.txt',delimiter='\t',names=['Names','Asymmetries'], overwrite = True)	


slope, intercept, r_value, p_value, std_err = stats.linregress(asymmetries,bovera)

print(r_value)

plt.scatter(asymmetries,bovera)

plt.show()

