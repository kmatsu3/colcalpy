#!/usr/bin/env python
# coding: UTF-8
import csv
import numpy
import sys
import math
import pylab

inputfile = csv.reader(open("parameter_input.txt", "r"), delimiter=":", quotechar='"')
for row in inputfile:
	if row[0] == "system_size":
		print(row[0],":",row[1])
		system_size = int(row[1])
	if row[0] == "tics":
		print(row[0],":",row[1])
		tics = int(row[1])
		if tics == 0:
			print("error 0 tics!")
			sys.exit()
bins=int(float(system_size)/float(tics))
histgram = numpy.zeros(bins)
correlation = numpy.zeros(bins)
spins = numpy.empty((0,2),float)
locations = numpy.empty((0,2),float)
count = 0
csv_reader = csv.reader(open("tmp.csv", "r"), delimiter=",", quotechar='"')
average_spin = numpy.zeros(2)
for row in csv_reader:
	count=count+1
	float_data = numpy.array(list(map(float,row)))
	spin = numpy.array([float_data[2]-float_data[0],float_data[3]-float_data[1]])
	norm = numpy.linalg.norm(spin)
	spin = numpy.array([spin/norm])
	spins = numpy.append(spins, spin, axis=0)
	locations = numpy.append(locations, numpy.array([[float_data[0],float_data[1]]]), axis=0)
for cell_1 in range(0,count):
	average_spin = average_spin + spins[cell_1,:]
	for cell_2 in range(0,count):
		product = numpy.dot(spins[cell_1,:],spins[cell_2,:])
		length = int(round(numpy.linalg.norm((locations[cell_1,:]-locations[cell_2,:])/float(tics)),0))
		if length < histgram.size :
			histgram[length] = histgram[length]+1.0
			correlation[length] = correlation[length] +product
		shift=numpy.array([system_size,0.0])
		length = int(round(numpy.linalg.norm((locations[cell_1,:]-locations[cell_2,:]+shift)/float(tics)),0))
		if length < histgram.size :
			histgram[length] = histgram[length]+1.0
			correlation[length] = correlation[length] +product
		length = int(round(numpy.linalg.norm((locations[cell_1,:]-locations[cell_2,:]-shift)/float(tics)),0))
		if length < histgram.size :
			histgram[length] = histgram[length]+1.0
			correlation[length] = correlation[length] +product
		shift=numpy.array([0.0,system_size])
		length = int(round(numpy.linalg.norm((locations[cell_1,:]-locations[cell_2,:]+shift)/float(tics)),0))
		if length < histgram.size :
			histgram[length] = histgram[length]+1.0
			correlation[length] = correlation[length] +product
		length = int(round(numpy.linalg.norm((locations[cell_1,:]-locations[cell_2,:]-shift)/float(tics)),0))
		if length < histgram.size :
			histgram[length] = histgram[length]+1.0
			correlation[length] = correlation[length] +product
		shift=numpy.array([system_size,system_size])
		length = int(round(numpy.linalg.norm((locations[cell_1,:]-locations[cell_2,:]+shift)/float(tics)),0))
		if length < histgram.size :
			histgram[length] = histgram[length]+1.0
			correlation[length] = correlation[length] +product
		length = int(round(numpy.linalg.norm((locations[cell_1,:]-locations[cell_2,:]-shift)/float(tics)),0))
		if length < histgram.size :
			histgram[length] = histgram[length]+1.0
			correlation[length] = correlation[length] +product
		shift=numpy.array([system_size,-system_size])
		length = int(round(numpy.linalg.norm((locations[cell_1,:]-locations[cell_2,:]+shift)/float(tics)),0))
		if length < histgram.size :
			histgram[length] = histgram[length]+1.0
			correlation[length] = correlation[length] +product
		length = int(round(numpy.linalg.norm((locations[cell_1,:]-locations[cell_2,:]-shift)/float(tics)),0))
		if length < histgram.size :
			histgram[length] = histgram[length]+1.0
			correlation[length] = correlation[length] +product
#		print(cell_1,cell_2,product,length)
normalization_factor=correlation[0]/histgram[0]-numpy.dot(average_spin,average_spin)/pow(float(count),2)
for bin in range(0,bins):
	if histgram[bin] > 0:
 		correlation[bin]=(correlation[bin]/histgram[bin]-numpy.dot(average_spin,average_spin)/pow(float(count),2))/normalization_factor
	if histgram[bin] <= 0:
 		correlation[bin]=0.0
#	print(float_data[0],float_data[1],float_data[2]-float_data[0],float_data[3]-float_data[1])
print(correlation)
print("origin", normalization_factor)
output = open("output_correlation.txt", "w")
for bin in range(0,bins):
	output.write(str(bin*tics) + " " + str(correlation[bin]/normalization_factor) + "\n")
output.close()
pylab.xlim([0,100/4])
pylab.xlabel("r/r_cell")
pylab.ylim([-1,1])
pylab.ylabel("g_s(r)/g_s(0)")
pylab.plot(correlation)
pylab.show()