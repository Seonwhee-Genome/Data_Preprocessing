import os, re, sys, stat
import json, pickle
import logging
import logging.handlers
import subprocess as sp

import pandas as pd
from GS_variants import GS_variant_effect_prediction
from glob import glob
from operator import methodcaller

def exec_DNA(date, symbol):
    fix = File_Rename(date, symbol)
    fix.rename_EQL2_fastq()
    fix.rename_EQL7_individual()
    fix.rename_EQL7_mutation()
    fix.rename_EQL2_copynumber()
    Fix = Directory_Rename(date, symbol)
    Fix.rename_EQL2_copynumber()
    Fix.rename_EQL7_mutation()
    Fix.rename_EQL7_individual()
    fix.link_from_EQL7_to_EQL2_fqgz()
    fix.link_from_EQL7_EQL2_to_EQL7_bam()

class Directory_Rename(object):
    
    def __init__(self, date, prefix):
        self.date = date
        self.prefix = prefix
        self.service_center = Rename_Center(date)

    def Dirs_belonging_to_Root_Dir(self, InputDir, OutputDir):
        DirNameList_1 = glob("%s*" % InputDir, OutputDir)
        DirNameList_2 = []
        for DirName in DirNameList_1:
            DirNameList_2.append(OutputDir + DirName.split("/")[-1])
        print("root directory is %s " % (OutputDir))
        self.Rename_process(DirNameList_2, OutputDir)
    
    def Rename_process(self, InputList, location):
        FirstFilter = list(map(methodcaller("split", "/"), InputList))  # e.g : [['', 'EQL2', 'GS_20170801', 'WXS', 'fastq', 'link', 'IRCR_BT16_1141_T_GS.fq.gz'], [.....], [......]]
        SecondFilter = list(map(lambda x: x[-1], FirstFilter))  # e.g : ['IRCR_BT16_1141_T_GS.fq.gz', '...']
        ThirdFilter = self.service_center.Sample_naming_inspection(SecondFilter)  # e.g : ['IRCR_BT16_1141_T_GS.fq.gz']
        ThirdFilter = list(set(ThirdFilter))
        CorretedFileNameList = self.service_center.heuristic_correction_of_naming_mistakes(ThirdFilter)
        print(CorretedFileNameList)
        OLDFILES = list(map(lambda x: location + x, ThirdFilter))
        NEWFILES = list(map(lambda x: location + x, CorretedFileNameList))
        if len(OLDFILES) == len(NEWFILES):
            for i in range(len(OLDFILES)):
                self.service_center.Rename(OLDFILES[i], NEWFILES[i])

    def rename_EQL7_individual(self):
        print("alphaA")
        DirNameList_before = glob("%s*" % self.service_center.bam_dir)
        print("root directory is %s " % (self.service_center.bam_dir))
        # FileNameList = glob("%s*_BS*" % (self.service_center.bam_dir))
        # FileNameList = glob("%s*_GS*" % (self.service_center.bam_dir))
        self.Rename_process(DirNameList_before, self.service_center.bam_dir)

    def rename_EQL7_mutation(self):
        print("alphaB")
        self.Dirs_belonging_to_Root_Dir(self.service_center.bam_dir, self.service_center.GS_test)

    def rename_EQL2_copynumber(self):
        print("alphaC")
        self.Dirs_belonging_to_Root_Dir(self.service_center.bam_dir, self.service_center.GBMscan)


class File_Rename(Directory_Rename):

    def __init__(self, date, prefix):
        super().__init__(date, prefix)

    def Files_belonging_to_A_Dir(self, InputDir, OutputDir):
        DirNameList_1 = glob("%s*" % InputDir)  # e.g : ['/EQL8/pipeline/SGI20170718/IRCR_BT16_1021_02_RSq/IRCR_BT16_1021_02_RSq_splice.bam', '...']
        DirNameList_2 = []
        for DirName in DirNameList_1:
            DirNameList_2.append(OutputDir + DirName.split("/")[-1])
        for oneDir in DirNameList_2:
            FileNameList = glob("%s/*" % oneDir)
            print("root directory is %s " % (oneDir))
            for oneFile in FileNameList:
                location = oneFile.split("/")[:-1]
                location = "/".join(location) + "/"
                self.Rename_process(FileNameList, location)

    def rename_EQL2_fastq(self):
        print("betaE")
        fastqs = glob("%s*" % self.service_center.fastq_link)
        self.Rename_process(fastqs, self.service_center.fastq_link)
        # e.g : [['', 'EQL2', 'GS_20170801', 'WXS', 'fastq', 'link', 'IRCR_BT16_1141_T_GS.fq.gz'], [.....], [......]]
        # e.g : ['IRCR_BT16_1141_T_GS.fq.gz', '...']
        # e.g : ['IRCR_BT16_1141_T_GS.fq.gz']
        
    def rename_EQL7_individual(self):
        print("gammaA")
        self.Files_belonging_to_A_Dir(self.service_center.bam_dir, self.service_center.bam_dir)

    def rename_EQL7_mutation(self):
        print("gammaB")
        self.Files_belonging_to_A_Dir(self.service_center.bam_dir, self.service_center.GS_test)

    def rename_EQL2_copynumber(self):
        print("gammaC")
        self.Files_belonging_to_A_Dir(self.service_center.bam_dir, self.service_center.GBMscan)

    def link_from_EQL7_to_EQL2_fqgz(self):
        print("betaA")
        self.service_center.Link_modifier('Fastq', '', 'GS')

    def link_from_EQL7_EQL2_to_EQL7_bam(self):
        print("betaC", "betaD")
        self.service_center.Link_modifier('BAM', '', 'GS')



class Rename_Center(object):

    def __init__(self, group, type="GS"):

        if type == "GS":
            self.groupName = "GS_" + group
            self.basedir = '/EQL7/pipeline/'
            self.bam_dir = '%s%s/' % (self.basedir, self.groupName)
            self.GS_test = "%sGS_test/" %(self.basedir)
            self.fastq_link = '/EQL2/%s/WXS/fastq/link/' %(self.groupName)
            self.GBMscan = '/EQL2/gbmscan_pool/'

        elif type == "BS":
            self.groupName = "BS_" + group
            self.basedir = '/EQL7/pipeline/'
            self.bam_dir = '%s%s/' % (self.basedir, self.groupName)
            self.BS_test = "%sBS_test/" % (self.basedir)
            self.fastq_link = '/EQL2/%s/WXS/fastq/link/' % (self.groupName)
            self.BTscan = '/EQL2/btscan_pool/'

        self.Pair = '_pair_filter_vep.dat'
        self.Single = '_single_filter_vep.dat'
        self.type = type

        Group_name = re.compile("^(GS|BS|SGI|MACROGEN|TERAGEN)_([0-9]{8})")
        if Group_name.match(self.groupName) != None:
            print("initialization done!")
        else:
            sys.exit(0)

    def Sample_naming_inspection(self, matchingList):
        standards = []
        needToRename = []

        if self.type == "GS":
            standard1 = re.compile("^S([0-9]{2})([0-9]{5,})(_T)_(GS)+")
            standard2 = re.compile("^NS([0-9]{2})_([0-9]{3,4})(|_T0[1-9])(|_P[0-9]{1,2}|_S[0-9]|_SP[0-9]{1,2})_(GS)+")
            standard3 = re.compile("^NS([0-9]{2})_([0-9]{3,4})(|_T0[1-9])(_BR[1-9]|_BR[1-9]S[0-9])(|_P[0-9]{1,2}|_SP[0-9]{1,2})_(GS)+")
            standard4 = re.compile("^(|HGF_)(IRCR|SNU|CHA|NCC|AMC)_([A-Z]{2,3})([0-9]{2})_([0-9]{3,4})(_B|_T|_T0[1-9])_(GS)+")
            standard5 = re.compile("^(|HGF_)(IRCR|SNU|CHA|NCC|AMC)_([A-Z]{2,3})([0-9]{2})_([0-9]{3,4})(|_T|_T0[1-9])(|_P[0-9]{1,2}|_S[0-9]|_SP[0-9]{1,2})_(GS)+")
            standard6 = re.compile("^(|HGF_)(IRCR|SNU|CHA|NCC|AMC)_([A-Z]{2,3})([0-9]{2})_([0-9]{3,4})(|_T0[1-9])(_BR[1-9]|_BR[1-9]S[0-9])(|_P[0-9]{1,2}|_SP[0-9]{1,2})_(GS|mm10_GS)+")
            standards = [standard1, standard2, standard3, standard4, standard5, standard6]
        elif self.type == "BS":
            standard1 = re.compile("^S([0-9]{2})([0-9]{5,})(_T)_(BS)+")
            standard2 = re.compile("^NS([0-9]{2})_([0-9]{3,4})(|_T0[1-9])(|_P[0-9]{1,2}|_S[0-9]|_SP[0-9]{1,2})_(BS)+")
            standard3 = re.compile("^NS([0-9]{2})_([0-9]{3,4})(|_T0[1-9])(_BR[1-9]|_BR[1-9]S[0-9])(|_P[0-9]{1,2}|_SP[0-9]{1,2})_(BS)+")
            standard4 = re.compile("^(|HGF_)(IRCR|SNU|CHA|NCC|AMC)_([A-Z]{2,3})([0-9]{2})_([0-9]{3,4})(_B|_T|_T0[1-9])_(BS)+")
            standard5 = re.compile("^(|HGF_)(IRCR|SNU|CHA|NCC|AMC)_([A-Z]{2,3})([0-9]{2})_([0-9]{3,4})(|_T|_T0[1-9])(|_P[0-9]{1,2}|_S[0-9]|_SP[0-9]{1,2})_(BS)+")
            standard6 = re.compile("^(|HGF_)(IRCR|SNU|CHA|NCC|AMC)_([A-Z]{2,3})([0-9]{2})_([0-9]{3,4})(|_T0[1-9])(_BR[1-9]|_BR[1-9]S[0-9])(|_P[0-9]{1,2}|_SP[0-9]{1,2})_(BS|mm10_BS)+")
            standards = [standard1, standard2, standard3, standard4, standard5, standard6]

        elif self.type == "WXS":
            standard1 = re.compile("^S([0-9]{2})([0-9]{5,})_T_(SS|TS)+")
            standard2 = re.compile("^NS([0-9]{2})_([0-9]{3,4})(|_T0[1-9])(|_P[0-9]{1,2}|_S[0-9]|_SP[0-9]{1,2})_(SS|TS)+")
            standard3 = re.compile("^NS([0-9]{2})_([0-9]{3,4})(|_T0[1-9])(_BR[1-9]|_BR[1-9]S[0-9])(|_P[0-9]{1,2}|_SP[0-9]{1,2})_(SS|TS)+")
            standard4 = re.compile("^(|HGF_)(IRCR|SNU|CHA|NCC|AMC)_([A-Z]{2,3})([0-9]{2})_([0-9]{3,4})(_B|_T|_T0[1-9])_(SS|TS)+")
            standard5 = re.compile("^(|HGF_)(IRCR|SNU|CHA|NCC|AMC)_([A-Z]{2,3})([0-9]{2})_([0-9]{3,4})(|_T0[1-9])(|_P[0-9]{1,2}|_S[0-9]|_SP[0-9]{1,2})_(SS|TS)+")
            standard6 = re.compile("^(|HGF_)(IRCR|SNU|CHA|NCC|AMC)_([A-Z]{2,3})([0-9]{2})_([0-9]{3,4})(|_T0[1-9])(_BR[1-9]|_BR[1-9]S[0-9])(|_P[0-9]{1,2}|_SP[0-9]{1,2})_(SS|TS|mm10_TS|mm10_SS)+")
            standards = [standard1, standard2, standard3, standard4, standard5, standard6]

        for ref in standards:
            needToRename = list(filter(lambda x: ref.match(x) == None, matchingList)) + needToRename

        return needToRename

    def heuristic_correction_of_naming_mistakes(self, inputList):
        outputList = []
        mistake1 = re.compile("^IRCR_([A-Z]{2,3})([0-9]{2})_([0-9]{3,4})(_0[1-9]|_0[1-9]_P0)_([A-Z]{2,3})+")
        mistake2 = re.compile("^S([0-9]{2})([0-9]{5,})_(GS)+")
        mistake3 = re.compile("^(|IRCR_)S([0-9]{2})_([0-9]{5,})_(GS)+")
        mistake4 = re.compile("^(|HGF_)(IRCR|SNU|CHA|NCC|AMC)_([A-Z]{2,3})([0-9]{2})_([0-9]{3,4})(_T)(_BR[1-9]|_BR[1-9]S[0-9])(|_P[0-9]{1,2}|_SP[0-9]{1,2})_(GS|mm10_GS)+")
        mistake5 = re.compile("^(RCR)_([A-Z]{2,3})([0-9]{2})_([0-9]{3,4})(_B|_T|_T0[1-9])_(GS)+")

        for NAME in inputList:
            if mistake1.match(NAME) != None:
                print("%s is Type 1 mistake; therefore _0 --> _T0" %(NAME))
                newNAME = NAME.replace("_0", "_T0")
            elif mistake2.match(NAME) != None:
                print("%s is Type 2 mistake; therefore _GS --> _T_GS" %(NAME))
                newNAME = NAME.replace("_GS", "_T_GS")
            elif mistake3.match(NAME) != None:
                if "IRCR_S" in NAME:
                    print("%s is Type 3 mistake; therefore IRCR_S --> S" %(NAME))
                    newNAME = NAME.replace("IRCR_S", "S")
                Delimiter = "_GS"
                newNAMEList = newNAME.split(Delimiter)
                newNAMEList[0] = ''.join(newNAMEList[0].split("_"))
                newNAME = Delimiter.join(newNAMEList)
                if mistake2.match(newNAME) != None:
                    print("%s is Type 2 mistake; therefore _GS --> _T_GS" %(NAME))
                    newNAME = newNAME.replace("_GS", "_T_GS")
            elif mistake4.match(NAME) != None:
                print("%s is Type 4 mistake; therefore _T --> " % (NAME))
                newNAME = NAME.replace("_T", "")
                print(newNAME)
            elif mistake5.match(NAME) != None:
                print("%s is Type 5 mistake; therefore RCR_ --> IRCR_" % (NAME))
                newNAME = NAME.replace("RCR_", "IRCR_")
                print(newNAME)
            else:
                print("No mistakes found : %s" % (NAME))
                newNAME = NAME
            outputList.append(newNAME)
        return outputList

    def Rename(self, A, B):
        os.system('mv %s %s' % (A, B))
        print("%s >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> \n %s " %(A, B))

    def Link_modifier(self, key, query="C00", type="RSq"):
        if type == "BS":
            Test = self.BS_test
        elif type == "GS":
            Test = self.GS_test
        Origin = {"Fastq": glob("%s*fq.gz" % (self.fastq_link)), "BAM": glob("%s*_%s/*recal*" % (self.bam_dir, type))}
        ToUnlink = {"Fastq": glob("%s*_%s/*fq.gz" % (self.bam_dir, type)),
                    "BAM": glob("%s*_%s/*recal*" % (Test, type)) + glob("%s*_%s/*recal*" % (self.GBMscan, type))}
        for unlink in ToUnlink[key]:
            target_matching = unlink.split("/")[-1]
            for origin in Origin[key]:
                target = origin.split("/")[-1]
                if target == target_matching:
                    try:
                        print("unlink ", unlink)
                        sp.check_call("unlink %s" % (unlink), shell=True)

                    except FileNotFoundError as fnf:
                        print("The symbolic link is already dead!!! %s" % (fnf))
                        sp.check_call("rm -rf %s" % (unlink), shell=True)
                        pass
                    print(origin, ">>>>>>>>>>>>>>", unlink)
                    os.symlink(origin, unlink)

    def FROM_GS_TO_BS_OR_VICE_VERSA(self, inputList):
        outputList = []
        for NAME in inputList:
            # newNAME = NAME.replace("_GS", "_BS")
            newNAME = NAME.replace("_BS", "_GS")
            outputList.append(newNAME)
            #sp.call("mv %s %s" % (NAME, newNAME), shell=True)
        return outputList
