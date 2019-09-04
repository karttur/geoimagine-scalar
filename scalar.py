'''
Created on 9 Nov 2018

@author: thomasgumbricht
'''

'''
Created on 9 Nov 2018

@author: thomasgumbricht
'''
#import os
#import numpy as np
from sys import exit
from geoimagine.ktnumba import ScalarTWIpercent
from geoimagine.mask import SingleBandMasking
import geoimagine.support.karttur_dt as mj_dt
#gc = garbage collect
import gc

class ProcessScalar:
    '''class for scalar processing'''   
    def __init__(self, process, session, verbose):
        self.session = session
        self.verbose = verbose
        self.process = process
        #print (self.process.proc.processid) 
  
        for locus in self.process.srcLayerD:
            for datum in self.process.srcLayerD[locus]:
                contFlag = True
                for dstcomp in self.process.dstLayerD[locus][datum]:
                    if self.process.dstLayerD[locus][datum][dstcomp]._Exists() and not self.process.overwrite:
                        self.session._InsertLayer(self.process.dstLayerD[locus][datum][dstcomp], self.process.overwrite, self.process.delete)
                        print ('scalar ready',self.process.dstLayerD[locus][datum][dstcomp].FPN)
                        continue 
                    for srccomp in self.process.srcLayerD[locus][datum]:
                        if not (self.process.srcLayerD[locus][datum][srccomp]):
                            print ('Src composition missing',datum)
                            contFlag = False
                            continue
                    if contFlag:
                        #Read the source layer
                        self.process.srcLayerD[locus][datum][srccomp].ReadRasterLayer()             
                        if self.process.proc.processid[0:10].lower() == 'twipercent': 
                            self._TWIpercent(locus,datum,srccomp,dstcomp)
                        elif self.process.proc.processid == 'convertdaytomonth':
                            self._ConvertDayToMonth(locus,datum,srccomp,dstcomp)
                        else:
                            exitstr = 'Unrecognized process in ProcessScalar %s' %(self.process.proc.processid)
                            exit(exitstr)
                         
                        #Postprocess the results   
                        SingleBandMasking(self.process.srcLayerD[locus][datum][srccomp], self.process.dstLayerD[locus][datum][dstcomp])
                        #copy the geoformat from the src layer
                        self.process.dstLayerD[locus][datum][dstcomp].CopyGeoformatFromSrcLayer(self.process.srcLayerD[locus][datum][srccomp].layer)
                        #write the results
                        self.process.dstLayerD[locus][datum][dstcomp].CreateDSWriteRasterArray()
                        #Register the layer
                        self.session._InsertLayer(self.process.dstLayerD[locus][datum][dstcomp], self.process.overwrite, self.process.delete)
                        #Set the class instances of the src and dst objects to None (using del causes the loop to crash)
                        self.process.dstLayerD[locus][datum][dstcomp].layer.NPBAND = None
                        self.process.srcLayerD[locus][datum][srccomp].layer.NPBAND = None
                        self.process.dstLayerD[locus][datum][dstcomp] = None
                        self.process.srcLayerD[locus][datum][srccomp] = None
                        #garbage collect to free memory
                        gc.collect()
                    
    def _TWIpercent(self,locus,datum,srccomp,dstcomp):
        '''
        '''
        self.process.dstLayerD[locus][datum][dstcomp].layer = lambda:None
        self.process.dstLayerD[locus][datum][dstcomp].layer.NPBAND = ScalarTWIpercent(
            self.process.srcLayerD[locus][datum][srccomp].layer.NPBAND, self.process.params.scalefac, 
            self.process.params.constant, self.process.params.divisor, self.process.params.power, 
            self.process.params.powfac, self.process.params.dstmax)

    def _ConvertDayToMonth(self, locus,datum,srccomp,dstcomp):
        dayD = {'01':31,'02':28,'03':31,'04':30,'05':31,'06':30,'07':31,'08':31,'09':30,'10':31,'11':30,'12':31}
        if len(datum) == 2:
            days = dayD[datum]
        elif len(datum) == 6:
            days = mj_dt.GetDaysInYYYY_MM(int(datum[0:4]), int(datum[4:6]))
        else:
            NOTYET
        self.process.dstLayerD[locus][datum][dstcomp].layer = lambda:None
        self.process.dstLayerD[locus][datum][dstcomp].layer.NPBAND = self.process.srcLayerD[locus][datum][srccomp].layer.NPBAND*self.process.params.factor*days + self.process.params.offset
 