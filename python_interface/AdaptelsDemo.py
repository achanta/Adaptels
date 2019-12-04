# This demo uses CFFI to complile and run the Adaptels algorithm

import os
import subprocess
from PIL import Image
from skimage.io import imread,imshow
import numpy as np
from timeit import default_timer as timer
from _adaptels.lib import Adaptels_main
from cffi import FFI



def segment(imgname,threshold,doRGBtoLAB):
	#--------------------------------------------------------------
	# read image and change image shape from (h,w,c) to (c,h,w)
	#--------------------------------------------------------------
	img = Image.open(imgname)
	# img = imread(imgname)
	img = np.asarray(img)
	print("input image shape", img.shape)

	dims = img.shape
	h,w,c = dims[0],dims[1],1
	if len(dims) > 1:
		c = dims[2]
		img = img.transpose(2,0,1)
	
	#--------------------------------------------------------------
	# Reshape image to a single dimensional vector
	#--------------------------------------------------------------
	img = img.reshape(-1).astype(np.double)
	labels = np.zeros((h,w), dtype = np.int32)
	numlabels = np.zeros(1,dtype = np.int32)
	#--------------------------------------------------------------
	# Prepare the pointers to pass to the C function
	#--------------------------------------------------------------
	ffibuilder = FFI()
	pinp = ffibuilder.cast("double*", ffibuilder.from_buffer(img))
	plabels = ffibuilder.cast("int*", ffibuilder.from_buffer(labels.reshape(-1)))
	pnumlabels = ffibuilder.cast("int*", ffibuilder.from_buffer(numlabels))

	
	start = timer()
	Adaptels_main(pinp,w,h,c,threshold,doRGBtoLAB,plabels,pnumlabels)
	end = timer()

	#--------------------------------------------------------------
	# Collect labels
	#--------------------------------------------------------------
	print("number of superpixels: ", numlabels[0])
	print("time taken in seconds: ", end-start)

	return labels.reshape(h,w),numlabels[0]


def drawBoundaries(imgname,labels,numlabels):

	img = Image.open(imgname)
	img = np.array(img)

	ht,wd = labels.shape

	for y in range(1,ht-1):
		for x in range(1,wd-1):
			if labels[y,x-1] != labels[y,x+1] or labels[y-1,x] != labels[y+1,x]:
				img[y,x,:] = 0

	return img

#------------------------------------------------------------------
# Before calling this function, compile the C code using
# "python compile_adaptels_lib.py"
#------------------------------------------------------------------
def adaptelsdemo():
	#--------------------------------------------------------------
	# Set threshold
	#--------------------------------------------------------------
	T = 60.0 # This value is best suited for 3-channel sRGB images
	doRGBtoLAB = True # only works if it is a three channel image
	imgname = "bee.png"

	labels,numlabels = segment(imgname,T,doRGBtoLAB)
	#--------------------------------------------------------------
	# Display segmentation result
	#------------------------------------------------------------
	segimg = drawBoundaries(imgname,labels,numlabels)
	# Image.fromarray(segimg).show()
	Image.fromarray(segimg).save("bee_adaptels.png")
	return

adaptelsdemo()



