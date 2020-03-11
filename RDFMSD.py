#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 10 12:28:24 2020

@author: s0199195 
"""

import matplot.lib.pyplot as plot #for plotting
import numpy as np #for maths
import pandas as pd #for numerical analysis, check how to use and use if it makes life easier.
import particle3D #for particle objects

#I'm going to have to calculate bulk density for the liquid state.  
#This is always the same for argon or any material. g(r)=(dnr/(4pidr*p)). 
#dnr is a function which computes numnber of particles in a shell thickness of dr.
#will have to figure out shell thickness
#for argon liquid there will be 12 neighbours to every molecule. coordination number 12, see wiki.

#for the solid part

def Binning(VMDposition, bin, box):
    #Reads in the array from VMD output file of all particle positions for each timestep and 
    #calculates the distances between particles and bins these into suitable bins
    #Normalise the results by number of particles.
    
    Nm = len(VMDposition) #defines the length between two particles 
    diff_array = [] #creates an array of the difference
    for i in range(Nm):
        for j in range(i):
            #finds distances between all pairs of particles 
            #and appends array.
            diff = VMDposition(j)-VMDposition(i) 
            cube1 = np.mod(diff,box)
            MICvecdiff = np.mod(cube1+(box/2),box)-(box/2)
            diff_array.append(np.linalg.norm(MICvecdiff))
            
    binValues = np.histogram(diff_array, bins=bins)[0]/(N/2)
    return binValues

def RDF(VMDposition, start, end, bins, box):
    #Takes an array of particle positions per each timestep and calculates an 
    #average from the start to end of simulation using the binning program above.
    
    Output = [BinnedParticles(pos, bins, box) for pos in VMDposition[start:end]

    RDFHistogram = sum(Output)/len(Output)
    volume = 4*pi*((bins[:-1]+bins[1:])/2)**2*(bins[1]-bins[0])
    RDFHistogram = RDFHistogram/volume
    RDFpoints = (bins[:-1]+bins[1:])/2
    
    return RDFHistogram, RDFpoints

def MSD(VMDpositions, start, end, box):
    
    inputted = VMDposition[start]
    MeanSD = []
    for time in range(start,end+1):
        tint = VMDposition[time]
        sum = 0
        for i in range(len(inputted)):
            diff = tint-inputted
            cube1 = np.mod(diff,box)#image in first cube
            MICvecdiff = np.mod(cube1+box/2,box)-box/2
            sum +=np.linalg.norm(MICvecdiff)**2
            
        MeanSD.append(sum/pos.shape[1])
    
    return MeanSD