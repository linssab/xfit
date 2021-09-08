import numpy as np
import os, sys
from .SimpleMath import *
from .FitMath import fit_single_spectrum

class Spectrum():
    def __init__(__self__,array=None,file_path=None):
        __self__.specfile = file_path
        if file_path.endswith(".mca"): 
            __self__.CALTAG = "<<CALIBRATION>>"
            __self__.DATTAG = "<<DATA>>"
        elif file_path.endswith(".spt"):
            __self__.CALTAG = "[CALIBRATION]"
            __self__.DATTAG = "[DATA]"
        if array is not None:
            if isinstance(array,np.ndarray):
                __self__.data = array.astype(np.float32)
            else:
                raise TypeError("Array must be a numpy.ndarray object!")
                return
        else:
            __self__.getdata()
        __self__.continuum = np.zeros(__self__.data.shape[0],dtype=np.float32)
        __self__.energyaxis = np.zeros(__self__.data.shape[0],dtype=np.float32)
        __self__.config = float
        __self__.gain = float
        __self__.FN = []

    def getdata(__self__):
        """ Extract the data contained in spectrum files

        ------------------------------------------------

        INPUT:
            mca; path
        OUTPUT:
            Data; 1D-array """
        
        mca = __self__.specfile
        name = str(mca)
        name = name.split("\\")[-1]
        name = name.replace('_',' ')

        # custom MC generated files
        if '#XRMC#' in name:
            Data = []
            datafile = open(mca, "r")
            lines = datafile.readlines()
            for line in lines:
                line = line.split()
                try: counts = float(line[1])
                except: counts = float(line[0])
                counts = counts * 10e3
                Data.append(counts)
            Data = np.asarray(Data)

        # this works for mca extension files
        else:
            datalist=[]
            datafile = open(mca, "r")
            line = datafile.readline()
            line = line.replace("\r","")
            line = line.replace("\n","")

            # AMPTEK files start with this tag
            if "<<PMCA SPECTRUM>>" in line:
                while "<<DATA>>" not in line:
                    line = datafile.readline()
                    if line == "": break
                line = datafile.readline()
                while "<<END>>" not in line:
                    try: datalist.append(int(line))
                    except ValueError as exception:
                        datafile.close()
                        raise exception.__class__.__name__
                    line = datafile.readline()
                    if line == "": break

            # Works if file is just counts per line
            elif line.isdigit():
                while line:
                    datalist.append(int(line))
                    line = datafile.readline()
                    if line == "": break
            
            elif "[LIDAQ1ch Spectrum]" in line:
                while "[DATA]" not in line:
                    line = datafile.readline()
                    if line == "": break
                line = datafile.readline()
                while line:
                    datalist.append(int(line))
                    line = datafile.readline()
                    if line == "": break

            # if file has two columns separated by space or tab
            elif "\t" in line or " " in line:
                while line:
                    counts = line.split("\t")[-1]
                    if counts.isdigit(): datalist.append(int(counts))
                    else:
                        counts = line.split(" ")[-1]
                        datalist.append(float(counts)*10e3)
                    line = datafile.readline()
                    if line == "": break
            datafile.close()
            Data = np.asarray(datalist, dtype=np.float32)
        __self__.data = Data
        return

    def get_anchors(__self__):
        mca_file = open(__self__.specfile,"r")
        line = mca_file.readline()
        param = []
        while line != "":
            while __self__.CALTAG not in line:
                line = mca_file.readline()
                if line == "": break
            while __self__.DATTAG not in line:
                line = mca_file.readline()
                if line == "": break
                line=line.replace('\r','')
                line=line.replace('\n','')
                line=line.replace('\t',' ')
                aux = line.split()
                try: param.append([int(aux[0]),float(aux[1])])
                except: pass 
            line = mca_file.readline()

        mca_file.close()
        if param == []: 
            raise ValueError("Couldn't fetch calibration from source!")
        for parameter in param:
            if len(param) <= 1: 
                raise ValueError("Need at least two calibration points!")
            elif parameter[0] < 0 or parameter[1] < 0: 
                raise ValueError("Cant receive negative values!")
            elif parameter[0] == 0 or parameter[1] == 0:
                raise ValueError("Cant receive zero as parameter!")
        else: pass
        return param
    
    def calibrate(__self__,x=None,y=None):
        if x is None and y is None:
            if __self__.specfile is None:
                raise Exception("No file was used as input!")
                return
            else: 
                params = __self__.get_anchors()
                if params == []: return
                x=[]
                y=[]
                for i in range(len(params)):
                    for k in range(len(params[i])):
                        x.append(params[i][0])
                        y.append(params[i][1])
        elif len(x) != len(y): 
            raise ValueError(f"X and Y must have the same dimension! Got {len(x)} and {len(y)}")
       
        __self__.gain, __self__.slope,  aa = list(linregress(x,y))
        for i in range(__self__.energyaxis.shape[0]):
            __self__.energyaxis[i] = __self__.gain*i + __self__.slope
    
    def fit_fano_and_noise(__self__):
        __self__.FN = FN_fit_gaus(
                __self__.data,
                __self__.continuum,
                __self__.energyaxis,
                __self__.gain)

    def create_pool(__self__,pool_file):
        __self__.pool = {}
        __self__.pool["bg"] = True
        __self__.pool["elements"] = {}
        f = open(pool_file, "r")
        lines = f.readlines()
        for line in lines:
            line = line.replace("\n","").replace("\r","").split("\t")
            element = line[0]
            __self__.pool["elements"][element] = []
            for peak in line[1:-1]:
                __self__.pool["elements"][element].append(peak)
        f.close()

    def fit(__self__, **kwargs):
        if "pool" in kwargs:
            pool = kwargs["pool"]
        else: pool = __self__.pool
        __self__.areas, __self__.plots = fit_single_spectrum(__self__, pool=pool)
