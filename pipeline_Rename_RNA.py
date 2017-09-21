import os, re, sys, stat
import json, pickle
import logging
import logging.handlers
import subprocess as sp

import pandas as pd
from GS_variants import GS_variant_effect_prediction
from glob import glob
from operator import methodcaller

def exec_RNA(date, symbol):
    File = File_Rename(date, symbol)
    #File.rename_EQL2_fastq()
    #File.rename_EQL8_expr()
    #File.rename_EQL8_mut()
    #Dir = Directory_Rename(date, symbol)
    #Dir.rename_EQL8_mut()
    #File.rename_EQL8_isoform()
    File.link_from_EQL8_to_EQL2_gsnap()
    #Dir.rename_EQL8_isoform()
    File.link_from_EQL8_to_EQL2_fqgz()
    #Dir.rename_EQL8_expr()


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
        FirstFilter = list(map(methodcaller("split", "/"),
                               InputList))  # e.g : [['', 'EQL2', 'GS_20170801', 'WXS', 'fastq', 'link', 'IRCR_BT16_1141_T_GS.fq.gz'], [.....], [......]]
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

    def rename_EQL8_expr(self):
        print("alphaA")
        DirNameList_before = glob("%s*" % self.service_center.tdf_dir)
        print("root directory is %s " % (self.service_center.tdf_dir))
        self.Rename_process(DirNameList_before, self.service_center.tdf_dir)
        # e.g : [['', 'EQL8', 'pipeline', 'SGI20170718', 'IRCR_BT16_1021_02_RSq'], [...]]
        # e.g : ['IRCR_BT16_1021_02_RSq', '..']
        # e.g : ['IRCR_BT16_1021_02_RSq']

    def rename_EQL8_mut(self):
        print("alphaB")
        DirNameList_before = glob("%s*" % self.service_center.bam_dir)
        print("root directory is %s " % (self.service_center.bam_dir))
        self.Rename_process(DirNameList_before, self.service_center.bam_dir)

    def rename_EQL8_isoform(self):
        print("alphaC")
        DirNameList_before = glob("%s*" % self.service_center.ei_dir)
        print("root directory is %s" % (self.service_center.ei_dir))
        self.Rename_process(DirNameList_before, self.service_center.ei_dir)
        DirNameList_before = glob("%s*" % self.service_center.skip_dir)
        print("root directory is %s" % (self.service_center.skip_dir))
        self.Rename_process(DirNameList_before, self.service_center.skip_dir)
        DirNameList_before = glob("%s*" % self.service_center.fusion_dir)
        print("root directory is %s" % (self.service_center.fusion_dir))
        self.Rename_process(DirNameList_before, self.service_center.fusion_dir)



class File_Rename(Directory_Rename):

    def __init__(self, date, prefix):

        Directory_Rename.__init__(self, date, prefix)

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
        print("betaA")
        fastqs = glob("%s*" % self.service_center.fastq_link)
        self.Rename_process(fastqs, self.service_center.fastq_link)
        # e.g : [['', 'EQL2', 'GS_20170801', 'WXS', 'fastq', 'link', 'IRCR_BT16_1141_T_GS.fq.gz'], [.....], [......]]
        # e.g : ['IRCR_BT16_1141_T_GS.fq.gz', '...']
        # e.g : ['IRCR_BT16_1141_T_GS.fq.gz']

    def rename_EQL8_expr(self):
        print("gammaA")
        self.Files_belonging_to_A_Dir(self.service_center.tdf_dir, self.service_center.tdf_dir)
        # e.g : ['/EQL8/pipeline/SGI20170718/IRCR_BT16_1021_02_RSq/IRCR_BT16_1021_02_RSq_splice.bam', '...']
        # e.g : [['', 'EQL8', 'pipeline', 'SGI20170718', 'IRCR_BT16_1021_02_RSq', 'IRCR_BT16_1021_02_RSq_splice.bam']]
        # e.g : ['IRCR_BT16_1021_02_RSq_splice.bam', 'IRCR_BT16_1141_RSq_splice.dedup.bai']
        # e.g : ['IRCR_BT16_1021_02_RSq_splice.bam']

    def rename_EQL8_mut(self):
        print("gammaB")
        self.Files_belonging_to_A_Dir(self.service_center.bam_dir, self.service_center.bam_dir)

    def rename_EQL8_isoform(self):
        print("gammaC")
        self.Files_belonging_to_A_Dir(self.service_center.ei_dir, self.service_center.ei_dir)
        self.Files_belonging_to_A_Dir(self.service_center.skip_dir, self.service_center.skip_dir)
        self.Files_belonging_to_A_Dir(self.service_center.fusion_dir, self.service_center.fusion_dir)

    def link_from_EQL8_to_EQL2_fqgz(self):
        print("betaB, betaC")
        self.service_center.Link_modifier('Fastq', '', 'RSq')

    def link_from_EQL8_to_EQL2_gsnap(self):
        print("betaD")
        self.service_center.Link_modifier('GSNAP', '', 'RSq')


class Rename_Center(object):
    def __init__(self, group):
        prefix = "SGI"
        self.groupName = prefix + "_" + group
        self.basedir = '/EQL8/pipeline/'
        self.all_dir = '%s%s%s_rsq2*/' % (self.basedir, prefix, group)
        self.tdf_dir = '%s%s%s_rsq2expr/' % (self.basedir, prefix, group)
        self.bam_dir = '%s%s%s_rsq2mut/' % (self.basedir, prefix, group)
        self.skip_dir = '%s%s%s_rsq2skip/' % (self.basedir, prefix, group)
        self.ei_dir = '%s%s%s_rsq2eiJunc/' % (self.basedir, prefix, group)
        self.fusion_dir = '%s%s%s_rsq2fusion/' % (self.basedir, prefix, group)
        self.fastq_link = '/EQL2/%s/RNASeq/fastq/link/' % (self.groupName)
        self.Pair = '_pair_filter_vep.dat'
        self.Single = '_single_filter_vep.dat'
        self.type = "RSq"

        Group_name = re.compile("^(SGI|MACROGEN|TERAGEN)_([0-9]{8})")
        if Group_name.match(self.groupName) != None:
            print("initialization done!")
        else:
            sys.exit(0)

    def Sample_naming_inspection(self, matchingList):

        needToRename = []
        standard1 = re.compile("^IRCR_([A-Z]{2,3})([0-9]{2})_([0-9]{3,4})(|_T0[1-9]|_T0[1-9]_P0)_(RSq)+")
        standard2 = re.compile("^IRCR_([A-Z]{2,3})([0-9]{2})_([0-9]{3,4})(_C0[0-9]|_T0[1-9]_C0[0-9])_(RSq)+")  # Total cells, Tumor Cells
        standard3 = re.compile("^IRCR_([A-Z]{2,3})([0-9]{2})_([0-9]{3,4})(_TAM|_T0[1-9]_TAM)_(RSq)+")  # TAM
        standards = [standard1, standard2, standard3]

        for ref in standards:
            needToRename = list(filter(lambda x: ref.match(x) == None, matchingList)) + needToRename

        return needToRename
    #########################################################################################################
    def TAM_sample_naming_correction(self, inputList):
        outputList = []
        mistake_TAM = re.compile("^IRCR_([A-Z]{2,3})([0-9]{2})_([0-9]{3,4})(|_T0[1-9]|_T0[1-9]_P0)_(RSq)+")

        for NAME in inputList:
            if mistake_TAM.match(NAME) != None:
                newNAME = NAME.replace("_RSq", "_C00_RSq")
            else:
                newNAME = NAME
            outputList.append(newNAME)
        return outputList
    #########################################################################################################

    def heuristic_correction_of_naming_mistakes(self, inputList):
        outputList = []
        mistake1 = re.compile("^IRCR_([A-Z]{2,3})([0-9]{2})_([0-9]{3,4})(_0[1-9]|_0[1-9]_P0)_([A-Z]{2,3})+")
        mistake6 = re.compile("^IRCR_([A-Z]{2,3})([0-9]{2})_([0-9]{3,4})(_T|_T_P0)_(RSq)+")

        for NAME in inputList:
            if mistake1.match(NAME) != None:
                print("%s is Type 1 mistake; therefore _0 --> _T0" % (NAME))
                newNAME = NAME.replace("_0", "_T0")
            elif mistake6.match(NAME) != None:
                print("%s is Type 6 mistake; therefore _T --> " % (NAME))
                newNAME = NAME.replace("_T", "")
            else:
                print("No mistakes found : %s" % (NAME))
                newNAME = NAME
            outputList.append(newNAME)
        return outputList

    def Rename(self, A, B):
        os.system('mv %s %s' % (A, B))
        print("%s >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> \n %s " % (A, B))


    def Link_modifier(self, key, query="C00", type="RSq"):
        Origin = {"Fastq": glob("%s/*%s*fq.gz" % (self.fastq_link, query)),
                  "GSNAP": glob("%s*%s_%s/*splice.gsnap.gz" % (self.bam_dir, query, type))}
        ToUnlink = {"Fastq": glob("%s*%s/*fq.gz" % (self.bam_dir, query)) + glob("%s*%s/*fq.gz" % (self.tdf_dir, query)),
                    "GSNAP": glob("%s*%s_%s/*splice.gsnap.gz" % (self.skip_dir, query, type)) +
                             glob("%s*%s_%s/*splice.gsnap.gz" % (self.ei_dir, query, type)) +
                             glob("%s*%s_%s/*splice.gsnap.gz" % (self.fusion_dir, query, type))}
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

    def Rewrite_Textfiles(self, oldExp, newExp, LogFileList):

        for afile in LogFileList:
            lines = []
            with open(afile, 'r') as of:
                for line in of:
                    if oldExp in line:
                        print(line)
                        line = line.replace(oldExp, newExp)
                        print("REPLACE >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>")
                        print(line)
                    lines.append(line)
            with open(afile, 'w') as wf:
                for line in lines:
                    wf.write(line)
                print("%s File changed!" % afile)


    def Rewrite_Textfiles_by_pattern(self, caseNum):

        LogFileList = glob("%s*/*txt" % (self.all_dir))
        #LogFileList = glob("%s*349*html" % (self.all_dir))

        if caseNum == 1:
            wrongPattern = ["(IRCR_[A-Z]{2,3}[0-9]{2}_[0-9]{3,4}_T_RSq)", "_T"]
            rightPattern = ["(IRCR_[A-Z]{2,3}[0-9]{2}_[0-9]{3,4}(|_T0[1-9]|_T0[1-9]_P0)_RSq)", ""]

        for afile in LogFileList:
            os.system("cp %s %s.backup" %(afile, afile))
            lines = []
            with open(afile, 'r') as of:
                for line in of:
                    wrong_word = re.findall(wrongPattern[0], line)
                    if len(wrong_word[0]) >= 1:
                        print(wrong_word[0])
                        right_word = wrong_word[0].replace(wrongPattern[1], rightPattern[1])
                        check = re.fullmatch(rightPattern[0], right_word)
                        if check != None:
                            line = line.replace(wrong_word[0], right_word)
                            print("REPLACE >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>")
                            print(line)
                        lines.append(line)   #### for safety
                    else:
                        lines.append(line)  ##### Be careful its location !!!!!!

            with open(afile, 'w') as wf:
                for line in lines:
                    wf.write(line)
                print("%s File changed!" % afile)

    def Maintain_or_Restore(self, FileList, choice="M"):
        for afile in FileList:
            if afile.split(".")[-1] == "backup":
                if choice == "M":
                    os.system("rm -f %s")
                    print("Remove %s" %afile)
                elif choice == "R":
                    os.system("mv %s %s" % (afile, afile[:-7]))
                    print("Restore %s --> %s" % (afile, afile[:-7]))







if __name__ == "__main__":
    exec_RNA("20170601", "RSq")

    #trial = Rename_Center("20170607")
    #trial.Rewrite_Textfiles_by_pattern(1)

