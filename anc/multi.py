#!/usr/bin/env python
# encoding: utf-8
"""
Example that reads columns and plots gnuplot style
"""

########################## Import Modules ########################## 
import sys
import os
import pylab as plt
import csv
import re
########################## Parameters ########################## 

CSVDATAFILE = "zoomEwsb"
CSVPLOTFILENAME = "reason1.eps"

## Flags ##



########################## Main ########################## 

def main(): 
    # this method is better, uses csv module
    data = ReadCSV()
    PlotCSV(data)


def ReadCSV(): #gets data from comma separated value file
    csvdata = csv.DictReader(CSVFileCleaner(open(CSVDATAFILE,'r')),delimiter = ' ',skipinitialspace=True )
    return csvdata

def PlotCSV(data):
    label = "data"
  
    # set limits
    xmin = 0.0
    xmax = 6.0    
    ymin = 0.0
    ymax = 11.0

    # set point style
    ebargs = dict(mc='None',alpha=1,ms=2*2.5,
              capsize=1.25,elinewidth=.5,mew=.5)                 
    fw = 4.5                    # width
    fh = 4.5                    # height
    plt.rc('text',usetex=True)
    plt.rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
    plt.rc('figure',figsize=(fw,fh))
    plt.xlabel(r'$\tilde x$')   
    plt.ylabel(r'$y$ ')
    # do the plot
    x,y,dy = [],[],[]
    for d in data:
      print d
      x.append(float(d['x']))
      y.append(float(d['y']))
      dy.append(float(d['dy']))
    #print x,y,dy
    plt.errorbar(x=x,y=y,yerr=dy,fmt='x',c='b',mec='b',label=label,**ebargs)  
    plt.xlim(xmin,xmax)
    plt.ylim(ymin,ymax) 
    plt.legend(loc='upper left',shadow=True, fancybox=True,numpoints=1)
    plt.savefig(CSVPLOTFILENAME,bbox_inches='tight')
    plt.cla()


class CSVFileCleaner: # removes comments from csv file
    def __init__ (self, file_to_wrap, comment_starts=('#', '//', '-- ','!')):
        self.wrapped_file = file_to_wrap
        self.comment_starts = comment_starts
    def next (self):
        line = self.wrapped_file.next()
        while line.strip().startswith(self.comment_starts) or len(line.strip())==0:
            line = self.wrapped_file.next()
        return line
    def __iter__ (self):
        return self


if __name__ == '__main__':
	main()
