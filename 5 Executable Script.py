#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
from scipy.stats import gaussian_kde
import os
import re
import tkinter as tk
from tkinter import filedialog
from tkinter import messagebox

root = tk.Tk()
root.withdraw()


# In[23]:


def loadData(filePath):
    # Load the matrix, skipping the first 3 rows of text
    data = np.loadtxt(filePath, skiprows=3)

    # Flatten the 2D matrix into a 1D array
    dataFlat = data.flatten()

    # Remove nan values (Profilm sets nan to +/- 10^38 for some reason)
    dataFiltered = dataFlat[np.abs(dataFlat) < 1000000]

    return dataFiltered

def findPeaks(data, nBins = 100, showPlot=True):
    
    # Calculate a histogram from data
    counts, bin_edges = np.histogram(data, bins=nBins)
    binMidpoints = (bin_edges[:-1] + bin_edges[1:]) / 2

    # Create a gaussian curve of data to improve peak finding accuracy 
    kde = gaussian_kde(data, bw_method=0.05)
    bin_width = np.diff(bin_edges)[0]
    total_samples = len(data)
    y = kde(binMidpoints)*bin_width*total_samples

    # Find the peaks in the counts data
    peaks, _ = find_peaks(counts, distance=10, prominence=max(counts)*0.3)
    peaksY, _ = find_peaks(y, distance=10, prominence=max(y)*0.1)
    
    # Visualize peaks
    if (showPlot):
        plt.hist(data, bins=nBins, color='blue')
        plt.plot(binMidpoints, y, color='red', linewidth=2, label='Smoothed KDE')
        for pos in peaks:
            plt.axvline(binMidpoints[pos], color='blue', linestyle='--', label='Simple Peak')
        for pos in peaksY:
            plt.axvline(binMidpoints[pos], color='red', linestyle='--', label='Peak Position')
        plt.legend()
        plt.show()

    return binMidpoints[peaksY]

def readProfilmScanFile(filePath):
    #Stores the filenames
    filename = []

    #Stores the X and Y scan positions
    x = []
    y= []
    
    try:
        with open(filePath, 'r') as file:
            #Skip the first line
            next(file) 
            
            #Iterate through the remaining lines
            for line in file:
                tmp = line.strip().split(',')

                #Format the scan position to Profilm scan file name
                filename.append('(' + tmp[0] + ', ' + tmp[1] + ').txt')

                if tmp[0] not in x:
                    x.append(tmp[0])
                if tmp[1] not in y:
                    y.append(tmp[1])

            x = np.sort(np.array(x).astype(float))
            y = np.sort(np.array(y).astype(float))[::-1]

            #A table containing the scan positions and height
            z = np.zeros((y.size+1, x.size+1)).astype(float)
            z[0,0] = np.nan
            z[1:,0] = y
            z[0,1:] = x
        return filename, z

    except FileNotFoundError:
        print("The file was not found.")
        return []
    except Exception as e:
        print(f"An error occurred: {e}")
        return []

def analyzeProfilmUniformity():
    
    #Ask for Scan Positions File
    scanFilename = filedialog.askopenfilename(title="Select the Scan Positions file", 
                                          message="Select the Scan Positions file",
                                          filetypes=[("Text files", "*txt"), ("All files", "*.*")])
    
    #Ask for Profilm data directory
    dataDir = filedialog.askdirectory(title="Select the Folder Containing Your Scan Files",
                                     message="Select the Folder Containing Your Scan Files",
                                     initialdir=scanFilename)

    #Ask for save data file name
    saveFilename = filedialog.asksaveasfilename(title="Save your TSV results",
                                                message="Save your TSV results",
                                                defaultextension=".tsv",
                                                filetypes=[("TSV files", "*.tsv"), ("All files", "*.*")],
                                                initialdir=dataDir)
    
    filename, z = readProfilmScanFile(scanFilename)
    r = z[:,0]
    c = z[0,:]
    for i in filename:
        data = loadData(dataDir + '/' + i)
        peaks = findPeaks(data,100,0)
        height = peaks[1] - peaks[0]
        tmp = i.replace(".txt", "").replace("(","").replace(")","")
        tmp = tmp.split(",")
        x = np.float64(tmp[0])
        y = np.float64(tmp[1])
        yi = np.isclose(r,y, atol = 0.1)
        xi = np.isclose(c,x, atol = 0.1)
        z[yi,xi] = height

    np.savetxt(saveFilename, z, delimiter="\t")


# In[21]:


analyzeProfilmUniformity()

