import os, re, sys, stat
import json, pickle
import logging
import logging.handlers
import subprocess as sp

import pandas as pd
from GS_variants import GS_variant_effect_prediction
from glob import glob
from operator import methodcaller



class Senior_Inspector(object):
    def __init__(self, group, type="GS"):

        if type == "GS":
            self.groupName = "GS_" + group
            self.basedir = '/EQL7/pipeline/'
            self.bam_dir = '%s%s/' % (self.basedir, self.groupName)
            self.GS_test = "%sGS_test/" % (self.basedir)
            self.fastq_link = '/EQL2/%s/WXS/fastq/link/' % (self.groupName)
            self.GBMscan = '/EQL2/gbmscan_pool/'

        elif type == "BS":
            self.groupName = "BS_" + group
            self.basedir = '/EQL7/pipeline/'
            self.bam_dir = '%s%s/' % (self.basedir, self.groupName)
            self.BS_test = "%sBS_test/" % (self.basedir)
            self.fastq_link = '/EQL2/%s/WXS/fastq/link/' % (self.groupName)
            self.BTscan = '/EQL2/btscan_pool/'


        elif type == "RSq":
            self.groupName = "MACROGEN_" + group + "_rsq2"
            # self.groupName = "SGI" + group + "_rsq2"
            self.basedir = '/EQL8/pipeline/'
            self.tdf_dir = '%s%sexpr/' % (self.basedir, self.groupName)
            self.bam_dir = '%s%smut/' % (self.basedir, self.groupName)
            self.skip_dir = '%s%sskip/' % (self.basedir, self.groupName)
            self.ei_dir = '%s%seiJunc/' % (self.basedir, self.groupName)
            self.fusion_dir = '%s%sfusion/' % (self.basedir, self.groupName)
            self.fastq_link = '/EQL2/SGI_%s/RNASeq/fastq/link/' % (group)

        self.Pair = '_pair_filter_vep.dat'
        self.Single = '_single_filter_vep.dat'
        self.type = type

    def file_existance(self, fileLists, Dir):

        List_Done = []
        List_Not_Done = []

        for aa in range(len(fileLists)):
            # ReferenceDirs = glob(Dir+"*"+fileLists[aa]+"*" + self.Pair)
            ReferenceDirs = glob(Dir + "*" + fileLists[aa] + "*")
            if ReferenceDirs == []:
                if not fileLists[aa].split('_')[-1] == "B":
                    List_Not_Done.append(fileLists[aa])
            else:
                print(ReferenceDirs)
                List_Done.append(fileLists[aa])

        print('\n\n\n\n\n', List_Done, " exists!!!")
        print('\n\n\n\n\n', List_Not_Done, " not exists")

    def inspection_for_Copynumber(self):
        CN = GS_Copynumber()
        CN.initialize_setting(self.groupName)
        CN.Select_TODO_list()

    def Duplicated_Dirs(self, dirName):
        for dir in dirName:
            path = "/EQL7/pipeline/GS_2*/%s*" % (dir)
            IsDuplicate = []
            searching = glob(path)
            for result in searching:
                IsADIR = os.path.isdir(result)
                if IsADIR == True:
                    IsDuplicate.append(result)
                else:
                    continue

            if len(IsDuplicate) >= 2:
                print("There are multiple %s, please check whether they were made by mistake " % (dir))
                print(IsDuplicate)

    def ReLink(self, toUnlink, newLink, target=''):
        STDOUT = sp.check_output("ls -l %s" % (toUnlink), shell=True)
        if target == '':
            TargetPoint = str(STDOUT).split("->")[-1][1:-3]
        else:
            TargetPoint = target

        try:
            sp.check_call("unlink %s" % (toUnlink), shell=True)

        except FileNotFoundError as fnf:
            print("The symbolic link is already dead!!! %s" % (fnf))
            sp.check_call("rm -rf %s" % (toUnlink), shell=True)
            pass
        os.symlink(TargetPoint, newLink)

    def Rename(self, A, B):
        os.system('mv %s %s' % (A, B))
        print("%s >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> \n %s " % (A, B))

    def Manual_Renaming(self, target_Dir):
        ## self.GS_test >> self.bam_dir
        olddir = "%sIRCR_GBM14_682_BR1_P5_RSq" % (target_Dir)
        newdir = "%sIRCR_GBM15_682_BR1_P5_RSq" % (target_Dir)
        matchingFiles = glob("%s/*" % (olddir))
        FirstFilter = list(map(methodcaller("split", "/"), matchingFiles))
        print("First Filter\n")
        print(FirstFilter)
        SecondFilter = list(map(lambda x: x[-1], FirstFilter))
        print("Second Filter\n")
        print(SecondFilter)
        ThirdFilter = list(set(SecondFilter))
        CorretedFileNameList = []
        print("Third Filter\n")
        print(ThirdFilter)

        for elements in ThirdFilter:
            filtered = elements.replace("_GBM14", "_GBM15")
            CorretedFileNameList.append(filtered)

        OLDFILES = list(map(lambda x: olddir + '/' + x, ThirdFilter))
        NEWFILES = list(map(lambda x: olddir + '/' + x, CorretedFileNameList))
        if len(OLDFILES) == len(NEWFILES):
            for i in range(len(OLDFILES)):
                self.Rename(OLDFILES[i], NEWFILES[i])
        self.Rename(olddir, newdir)

    def read_DAT(self, isPAIR=True):
        VEP = GS_variant_effect_prediction()
        VEP.initialize_setting(self.groupName)
        VEP.initialize_Single_Pair(isPair=isPAIR)
        LIST = VEP.Select_TODO_list()
        INDEX = ['Sample_Name', 'CHR', 'POSITION', 'REF', 'ALT', 'n_nRef', 'n_nAlt', 't_nRef', 't_nAlt', 'Annot_type',
                 'canonical', 'SC']
        instant_msg = []

        for dirs in LIST:
            ID_mutect = dirs.split('/')[-1] + '.mutect'
            ID_indels = dirs.split('/')[-1] + '.indels'
            for ID in [ID_mutect, ID_indels]:
                if isPAIR:
                    file = dirs + '/' + ID + self.Pair
                else:
                    file = dirs + '/' + ID + self.Single
                try:
                    with open(file, 'r') as f:
                        for line in f:
                            (SName, CHR, POS, REF, ALT, n_nRef, n_nAlt, t_nRef, t_nAlt, AnnotGene, AnnotCH_dna,
                             AnnotCH_prot, Annot_type, canonical, refseq, SC) = line.split('\t')
                            Sample = [SName, CHR, POS, REF, ALT, n_nRef, n_nAlt, t_nRef, t_nAlt, Annot_type, canonical,
                                      SC]
                            Sample_info = pd.Series(data=Sample, index=INDEX)
                            for element in Sample:
                                if element == '':
                                    print(line, Sample_info)
                                    instant_msg.append(file)
                                    with open("/home/jsgene/JK1/NGS/exec_history/%s_missing_in_DAT.txt" % ID,
                                              'a') as wf:
                                        wf.write(line)
                                else:
                                    print(dirs, 'is OK')
                                    continue
                except IOError as err:
                    logging.basicConfig(filename="/home/jsgene/JK1/NGS/exec_history/%s_no_DAT_file.log" % ID,
                                        level=logging.DEBUG)
                    logging.debug(err)
                    continue
            print("something wrong with %s" % instant_msg)


# def sequence_Quality(self):

###########################################################################################################################################
###########################################################################################################################################




if __name__=="__main__":
    ## Sample Names
    # sampleNames = glob("%s/*" % (self.bam_dir))
    # sampleNames = glob("%s*_BS"%(self.BS_test))
    ##sampleNames = glob("%s*_GS" % (self.BTscan))
    #sampleNames = glob("/EQL7/pipeline/BS_20170719/*BS")

    # sampleNames = glob("%s*_GS/*" % (self.BTscan))
    # sampleNames = glob("%s/*" % (self.BTscan))
    # sampleNames = ["%sIRCR_S16_05210_GS" % (self.GS_test)]
    # sampleNames = glob("/EQL8/pipeline/TERAGEN_20160711_rsq2mut/*") + glob("/EQL8/pipeline/TERAGEN_20160801_rsq2mut/*") + glob("/EQL8/pipeline/TERAGEN_20160826_rsq2mut/*") + glob("/EQL8/pipeline/TERAGEN_20160906_rsq2mut/*") + glob("/EQL8/pipeline/TERAGEN_20160907_rsq2mut/*") + glob("/EQL8/pipeline/TERAGEN_20160912_rsq2mut/*")
    # sampleNames = glob("/EQL8/pipeline/TERAGEN_20160801_rsq2expr/*") + glob("/EQL8/pipeline/TERAGEN_20160826_rsq2expr/*1076*") + glob("/EQL8/pipeline/TERAGEN_20160906_rsq2expr/*") + glob("/EQL8/pipeline/TERAGEN_20160907_rsq2expr/*") + glob("/EQL8/pipeline/TERAGEN_20160912_rsq2expr/*")
    # sampleNames = ["%sIRCR_BT16_1290_T02_P0_GS" %(self.GS_test)]

    ##sampleNames = glob("/EQL7/pipeline/GS_2017*/*_BR*")
    #sampleNames = glob("/EQL7/pipeline/GS_test/*BR*")
    ##sampleNames = glob("/EQL7/pipeline/GS_test/*_BR*")
    ###sampleNames = glob("/EQL2/GS_*/WXS/fastq/link/*_BR*")
    #####sampleNames = glob("/EQL7/pipeline/GS_20161206/*")
    sampleNames = glob("/EQL7/pipeline/GS_test/*")

    print(sampleNames)


    ############################ Change Dir names and File names #######################################################
    ins = Inspector("GS_20161206", "GS")
    #ins = Inspector("20170830", "RSq")

    #LogFileList = glob("/EQL8/pipeline/TERAGEN_20160711_rsq2expr/*C00_RSq/*log") + glob("/EQL8/pipeline/TERAGEN_20160801_rsq2expr/*C00_RSq/*log") + glob("/EQL8/pipeline/TERAGEN_20160826_rsq2expr/*1076*C00_RSq/*log") + glob("/EQL8/pipeline/TERAGEN_20160906_rsq2expr/*C00_RSq/*log") + glob("/EQL8/pipeline/TERAGEN_20160907_rsq2expr/*C00_RSq/*log") + glob("/EQL8/pipeline/TERAGEN_20160912_rsq2expr/*C00_RSq/*log")
    LogFileList = glob("/EQL8/pipeline/TERAGEN_20160*_rsq2*/*html")
    # LogFileList = glob("/EQL7/pipeline/BS_*/*_BS/*log")
    #### LogFileList = glob("/EQL2/btscan_pool/BS_*/*_BS/*log") + glob("/EQL7/pipeline/BS_*/*_BS/*log")
    #LogFileList = glob("/EQL8/pipeline/SGI20170607*/*/*log")
    #LogFileList = glob("/EQL8/pipeline/SGI20170601*/*html")
    ins.read_Logs("IRCR_BT16_1081_T01", "IRCR_BT16_1081_T01_C00", LogFileList)
    ins.read_Logs("IRCR_BT16_1077", "IRCR_BT16_1077_C00", LogFileList)
    ins.read_Logs("IRCR_BT16_1085_T01", "IRCR_BT16_1085_T01_C00", LogFileList)
    ins.read_Logs("IRCR_BT16_1094_T01", "IRCR_BT16_1094_T01_C00", LogFileList)
    ins.read_Logs("IRCR_BT16_1088_T01", "IRCR_BT16_1088_T01_C00", LogFileList)
    ins.read_Logs("IRCR_BT16_1076_T01", "IRCR_BT16_1076_T01_C00", LogFileList)
    ins.read_Logs("IRCR_BT16_1076_T02", "IRCR_BT16_1076_T02_C00", LogFileList)
    ins.read_Logs("IRCR_BT16_1036_T03", "IRCR_BT16_1036_T03_C00", LogFileList)


    #ins.read_Logs("IRCR_LC14_342_T_RSq", "IRCR_LC14_342_RSq", LogFileList)  ## Editing Log Files
    #ins.read_Logs("IRCR_LC14_346_T_RSq", "IRCR_LC14_346_RSq", LogFileList)
    #ins.read_Logs("IRCR_LC14_370_T_RSq", "IRCR_LC14_370_RSq", LogFileList)
    #ins.read_Logs("IRCR_LC15_446_T_RSq", "IRCR_LC15_446_RSq", LogFileList)
    #ins.read_Logs("IRCR_LC14_337_T_RSq", "IRCR_LC14_337_RSq", LogFileList)
    #ins.read_Logs("IRCR_LC14_349_T_RSq", "IRCR_LC14_349_RSq", LogFileList)

    #ins.read_Logs("GS_", "BS_")  ## Editing Log Files
    ###ins.Manual_Renaming(ins.bam_dir)
    ##inputs = glob("/EQL7/pipeline/BS_2017*/*_GS")
    ##ins.change_GS_dir_to_BS_dir(inputs)

    ############################# Reset Symbolic Links and Create again ################################################
    # ins.ReLink('/EQL7/pipeline/GS_test/S16124067_T_GS/S16124067_T_GS.recal.bai', '/EQL7/pipeline/GS_test/S16124067_T_GS/S16124067_T_GS.recal.bai', '/EQL7/pipeline/GS_20170711/S16124067_T_GS/S16124067_T_GS.recal.bai')


    ########################### Find out duplicated dirs ###############################################################
    # ins.Duplicated_Dirs(["IRCR_BT16_1061_B", "IRCR_BT16_1061_T", "IRCR_BT16_943_B", "IRCR_BT16_943_T", "IRCR_BT15_842_B", "IRCR_BT15_842_T", "IRCR_BT16_1165_B", "IRCR_BT16_1165_T", "IRCR_GBM12_140_BR1", "IRCR_BT16_966_B", "IRCR_BT16_966_T02", "IRCR_BT16_966_T04", "IRCR_GBM15_717_B", "IRCR_GBM15_717_T", "IRCR_BT16_936_B", "IRCR_BT16_936_T01", "IRCR_BT16_936_T02", "IRCR_BT15_884_B", "IRCR_BT15_884_T01", "IRCR_BT15_884_T05", "IRCR_BT16_1047_B", "IRCR_BT16_1047_T", "IRCR_BT15_890_B", "IRCR_BT15_890_T", "IRCR_BT16_937_B", "IRCR_BT16_937_T", "IRCR_BT16_1045_B", "IRCR_BT16_1045_T01", "IRCR_BT16_1045_T02", "IRCR_BT16_1045_T03", "IRCR_BT16_1045_T05", "IRCR_BT16_1072_B", "IRCR_BT16_1072_T01", "IRCR_BT16_1072_T02", "IRCR_BT16_1081_B", "IRCR_BT16_1081_T02", "IRCR_BT16_1105_B", "IRCR_BT16_1105_T01", "IRCR_BT16_1105_T02", "IRCR_BT16_948_B", "IRCR_BT16_948_T01", "IRCR_BT16_948_T02", "IRCR_BT16_1112_B", "IRCR_BT16_1112_T", "IRCR_BT16_945_B", "IRCR_BT16_945_T", "IRCR_BT16_951_B", "IRCR_BT16_951_T01", "IRCR_BT16_951_T03", "IRCR_BT16_951_T05", "IRCR_BT16_1059_B", "IRCR_BT16_1059_T01", "IRCR_BT16_1059_T02", "IRCR_BT16_1078_B", "IRCR_BT16_1078_T", "IRCR_BT16_1141_B", "IRCR_BT16_1141_T", "IRCR_BT16_1144_B", "IRCR_BT16_1144_T", "IRCR_BT16_1146_B", "IRCR_BT16_1146_T", "IRCR_BT16_1147_B", "IRCR_BT16_1147_T", "IRCR_BT16_1149_B", "IRCR_BT16_1149_T", "IRCR_BT16_1153_B", "IRCR_BT16_1153_T01", "IRCR_BT16_1153_T02", "IRCR_BT16_1153_T03", "IRCR_BT16_1160_B", "IRCR_BT16_1160_T01", "IRCR_BT16_1160_T02", "IRCR_BT16_1170_B", "IRCR_BT16_1170_T", "IRCR_BT16_1173_B", "IRCR_BT16_1173_T01", "IRCR_BT16_1173_T02", "IRCR_BT16_1174_B", "IRCR_BT16_1174_T", "IRCR_BT16_1179_B", "IRCR_BT16_1179_T", "IRCR_BT16_1185_B", "IRCR_BT16_1185_T", "IRCR_BT16_1197_B", "IRCR_BT16_1197_T", "IRCR_BT16_1199_B", "IRCR_BT16_1199_T", "IRCR_BT16_995_B", "IRCR_BT16_995_T", "IRCR_BT16_1140_B", "IRCR_BT16_1140_T", "IRCR_BT16_1143_B", "IRCR_BT16_1143_T", "IRCR_BT16_1217_B", "IRCR_BT16_1217_T", "IRCR_BT16_1205_B", "IRCR_BT16_1205_T01", "IRCR_BT16_1184_B", "IRCR_BT16_1184_T", "IRCR_BT16_1195_B", "IRCR_BT16_1195_T", "IRCR_BT16_1193_B", "IRCR_BT16_1193_T", "IRCR_BT16_1186_B", "IRCR_BT16_1186_T01", "IRCR_BT16_1186_T02", "IRCR_BT16_1208_B", "IRCR_BT16_1208_T", "IRCR_BT15_904_B", "IRCR_BT15_904_T01", "IRCR_BT15_904_T04", "IRCR_BT15_904_T05", "IRCR_BT16_1218_B", "IRCR_BT16_1218_T", "IRCR_BT16_1231_B", "IRCR_BT16_1231_T", "IRCR_BT16_1232_B", "IRCR_BT16_1232_T01", "IRCR_BT16_1232_T04", "IRCR_BT16_1236_B", "IRCR_BT16_1236_T", "IRCR_BT16_1237_B", "IRCR_BT16_1237_T", "IRCR_BT16_1240_B", "IRCR_BT16_1240_T", "IRCR_BT16_1245_B", "IRCR_BT16_1245_T02", "IRCR_BT16_1249_B", "IRCR_BT16_1249_T01", "IRCR_BT16_1249_T02", "IRCR_BT16_1250_B", "IRCR_BT16_1250_T", "IRCR_BT17_1258_B", "IRCR_BT17_1258_T02", "NCC_BT17_001_B", "NCC_BT17_001_T03", "NCC_BT17_001_T02", "NCC_BT17_001_T01", "IRCR_BT17_1269_B", "IRCR_BT17_1269_T01", "IRCR_BT17_1269_T02", "IRCR_BT17_1271_B", "IRCR_BT17_1271_T", "IRCR_BT17_1279_B", "IRCR_BT17_1279_T", "IRCR_BT17_1261_B", "IRCR_BT17_1261_T", "IRCR_BT17_1262_B", "IRCR_BT17_1262_T", "IRCR_BT17_1264_B", "IRCR_BT17_1264_T", "IRCR_BT17_1267_B", "IRCR_BT17_1267_T", "IRCR_BT17_1277_B", "IRCR_BT17_1277_T", "IRCR_BT17_1287_B", "IRCR_BT17_1287_T", "IRCR_BT17_1290_B", "IRCR_BT17_1290_T04", "IRCR_BT17_1290_T05", "IRCR_BT17_1290_T06", "IRCR_BT17_1291_B", "IRCR_BT17_1291_T02", "IRCR_BT17_1291_T03", "IRCR_BT17_1297_B", "IRCR_BT17_1297_T", "IRCR_BT17_1299_B", "IRCR_BT17_1299_T", "IRCR_BT17_1301_B", "IRCR_BT17_1301_T", "IRCR_BT17_1305_B", "IRCR_BT17_1305_T02", "IRCR_BT17_1305_T03", "IRCR_BT17_1305_T04", "CHA_BT17_002_B", "CHA_BT17_002_T", "SNU_BT17_003_B", "SNU_BT17_003_T", "SNU_BT17_008_B", "SNU_BT17_008_T", "SNU_BT17_009_B", "SNU_BT17_009_T", "NCC_BT17_011_B", "NCC_BT17_011_T02", "AMC_BT17_005_B", "AMC_BT17_005_T01", "AMC_BT17_006_B", "AMC_BT17_006_T", "AMC_BT17_007_B", "AMC_BT17_007_T", "IRCR_BT17_1303_B", "IRCR_BT17_1303_T02", "IRCR_BT17_1303_T03", "IRCR_BT17_1306_B", "IRCR_BT17_1306_T", "IRCR_BT17_1307_B", "IRCR_BT17_1307_T", "IRCR_BT15_891_B", "IRCR_BT15_891_T01", "IRCR_S13_18157", "IRCR_S16_05210", "IRCR_BT17_1319_B", "SNU_BT17_020_B", "SNU_BT17_020_T", "IRCR_BT17_1324_B", "IRCR_BT17_1324_T", "IRCR_BT17_1326_B", "IRCR_BT17_1326_T", "IRCR_BT17_1330_B", "IRCR_BT17_1330_T01", "IRCR_BT17_1330_T02", "IRCR_BT17_1330_T03", "IRCR_BT17_1330_T06", "IRCR_BT17_1332_B", "IRCR_BT17_1332_T", "IRCR_BT15_870_B", "IRCR_BT15_870_T", "IRCR_BT16_1023_B", "IRCR_BT16_1023_T", "IRCR_GBM13_300_B", "IRCR_GBM13_300_T", "LMJ_VIRUS_002_T", "IRCR_BT17_1317_B", "IRCR_BT17_1317_T01", "IRCR_BT17_1331_B", "IRCR_BT17_1331_T04", "IRCR_BT17_1339_B", "IRCR_BT17_1339_T01", "IRCR_BT17_1339_T02", "IRCR_BT17_1339_T03", "IRCR_BT17_1339_T04", "CHA_BT17_016_B", "CHA_BT17_016_T01", "CHA_BT17_016_T02", "CHA_BT17_016_T03", "CHA_BT17_016_T04", "SNU_BT17_013_B", "SNU_BT17_013_T", "SNU_BT17_014_B", "SNU_BT17_014_T", "SNU_BT17_015_B", "SNU_BT17_015_T", "SNU_BT17_017_B", "SNU_BT17_017_T", "IRCR_BT17_1272_B", "IRCR_BT17_1272_T", "IRCR_GBM14_400_B", "IRCR_GBM14_400_T", "IRCR_GBM15_766_B", "IRCR_GBM15_766_T", "IRCR_GBM14_608_BR1S2", "IRCR_BT16_935_B", "IRCR_BT16_935_T", "IRCR_BT15_925_B", "IRCR_BT15_925_T", "IRCR_BT17_1323_B", "IRCR_BT17_1323_T01", "IRCR_BT17_1323_T02", "IRCR_BT17_1323_T04", "IRCR_BT17_1325_B", "IRCR_BT17_1325_T01", "IRCR_BT17_1325_T03", "IRCR_BT17_1332_B", "IRCR_BT17_1332_T", "IRCR_GBM10_038_B", "IRCR_GBM10_038_T", "IRCR_GBM13_231_B", "IRCR_GBM13_231_T", "IRCR_GBM13_292_B", "IRCR_GBM13_292_T", "IRCR_GBM13_319_B", "IRCR_GBM13_319_T", "IRCR_GBM15_711_B", "IRCR_GBM15_711_T", "IRCR_GBM14_458_B", "IRCR_GBM14_458_T", "IRCR_GBM14_487_B", "IRCR_GBM14_487_T", "IRCR_GBM14_494_B", "IRCR_GBM14_494_T", "IRCR_GBM14_497_B", "IRCR_GBM14_497_T", "IRCR_GBM14_508_B", "IRCR_GBM14_508_T", "IRCR_GBM13_316_BR1S1", "IRCR_GBM15_698_BR1S2", "IRCR_BT15_849_T02_BR1S1", "IRCR_BT16_930_T02_BR1S1", "IRCR_GBM15_677_BR1S2", "IRCR_BT15_802_T03_BR1S1", "IRCR_GBM14_542_T01_BR1S1", "IRCR_BT15_891_T01_BR1S2", "IRCR_GBM14_589_B", "IRCR_GBM14_589_T", "IRCR_GBM14_567_B", "IRCR_GBM14_567_T", "IRCR_GBM14_586_B", "IRCR_GBM14_586_T", "IRCR_GBM14_626_B", "IRCR_GBM14_626_T", "IRCR_BT16_1134_T", "IRCR_BT16_1154_T", "IRCR_BT15_1134_B", "IRCR_BT16_1154_B", "IRCR_BT17_1290_T01", "IRCR_BT17_1342_B", "IRCR_BT17_1342_T", "IRCR_BT17_1343_B", "IRCR_BT17_1343_T", "IRCR_BT17_1346_B", "IRCR_BT17_1346_T01", "IRCR_BT17_1350_B", "IRCR_BT17_1350_T", "IRCR_BT17_1354_B", "IRCR_BT17_1354_T01", "IRCR_BT17_1354_T02", "IRCR_BT17_1357_B", "IRCR_BT17_1357_T", "SNU_BT17_022_B", "SNU_BT17_022_T01", "SNU_BT17_022_T02", "AMC_BT17_021_B", "AMC_BT17_021_T01", "AMC_BT17_021_T02", "IRCR_BT17_1303_T01", "IRCR_BT17_1353_B", "IRCR_BT17_1353_T01", "IRCR_BT17_1353_T02", "IRCR_BT17_1353_T03", "CHA_BT17_023_B", "CHA_BT17_023_T03", "CHA_BT17_023_T04", "SNU_BT17_026_B", "SNU_BT17_026_T01", "SNU_BT17_027_B", "SNU_BT17_027_T", "IRCR_BT17_1379_B", "IRCR_BT17_1379_T", "IRCR_BT17_1384_B", "IRCR_BT17_1384_T", "IRCR_BT17_1351_B", "S16124067_T", "SNU_BT17_024_B", "SNU_BT17_024_T01", "SNU_BT17_025_B", "SNU_BT17_025_T01", "IRCR_BT17_1366_B", "IRCR_BT17_1366_T", "IRCR_BT17_1367_B", "IRCR_BT17_1367_T", "IRCR_BT17_1373_B", "IRCR_BT17_1373_T01", "IRCR_BT17_1373_T02", "IRCR_BT17_1373_T04", "IRCR_BT17_1375_B", "IRCR_BT17_1375_T", "IRCR_BT17_1340_B", "IRCR_BT17_1340_T", "IRCR_BT16_1173_T06_P0", "IRCR_BT16_1245_T01_P0", "IRCR_BT17_1247_B", "IRCR_BT16_1247_T", "IRCR_BT16_1290_T02_P0", "IRCR_BT16_1177_T"])

    ########################### Inspect NGS output files whether their forms and dimensions interrupted or not #########
    # ins.read_DAT()

