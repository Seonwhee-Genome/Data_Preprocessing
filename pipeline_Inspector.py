import os, re, sys, stat
import json, pickle
import logging
import logging.handlers
import subprocess as sp

import pandas as pd
from GS_variants import GS_variant_effect_prediction
from glob import glob
from operator import methodcaller
from pipeline_Inspection_office import Senior_Inspector

class Junior_Inspector(Senior_Inspector):
    def Duplicated_Dirs(self):
        singleList = glob("/EQL7/pipeline/GS_single/IRCR_*BR*")
        for si in singleList:
            query = si.split("/")[-1]
            newpath = "/EQL7/pipeline/GS_test/%s" % query
            Silhwa = os.path.isdir(newpath)
            if Silhwa == True:
                print(si)
                print(newpath)


########################################################################################################################

if __name__=="__main__":
    ins = Junior_Inspector("GS_20161206", "GS")
    ins.Duplicated_Dirs()
