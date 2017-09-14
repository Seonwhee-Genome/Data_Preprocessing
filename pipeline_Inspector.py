import os, re, sys, stat
import json, pickle
import logging
import logging.handlers
import subprocess as sp

import pandas as pd
from GS_variants import GS_variant_effect_prediction
from glob import glob
from operator import methodcaller


class Logger(object):

    def __init__(self):
        self.Save_file_path = '/home/jsgene/JK1/NGS/exec_history/'  # log files that indicate execution time
        self.JSON_path = '/data1/home/jsgene/JSONS/'
        self.PICKLE_path = '/data1/home/jsgene/PICKLES'
        self.Desktop = '/home/jsgene/JK1/NGS/'

    def GS_TIME(self):
        now = time.localtime()
        begin_time = '%04d-%02d-%02d %02d:%02d:%02d' % (now.tm_year, now.tm_mon, now.tm_mday, now.tm_hour, now.tm_min, now.tm_sec)
        return begin_time

    def GS_exec_log(self, filename, saving_data):
        # Save history of operation
        with open('%s%s' % (self.Save_file_path, filename), 'a') as f:
            f.write(saving_data)

    def GS_progress_log(self, log_filename, saving_data):

        normal_logger = logging.getLogger()
        normal_formatter = logging.Formatter('[%(process)d %(processName)s %(thread)d %(threadName)s|%(filename)s:%(lineno)s] %(asctime)s > %(message)s')

        # create file handler which logs even debug messages
        normal_fileHandler = logging.FileHandler('%s%s.log' % (self.Save_file_path, log_filename))
        normal_streamHandler = logging.StreamHandler()

        # create formatter and add it to the handlers
        normal_streamHandler.setFormatter(normal_formatter)
        normal_fileHandler.setFormatter(normal_formatter)

        normal_logger.addHandler(normal_fileHandler)
        normal_logger.addHandler(normal_streamHandler)
        normal_logger.setLevel(logging.DEBUG)
        normal_logger.debug(saving_data)


    def GS_error_log(self, log_filename, function_name, message, exc_info=True):
        try:
            # execution of the function
            function_name
        except StandardError as e:
            stderr_logger = logging.getLogger()
            stderr_formatter = logging.Formatter('[%(process)d %(processName)s %(thread)d %(threadName)s|%(filename)s:%(lineno)s] %(asctime)s > %(message)s')

            stderr_fileHandler = logging.FileHandler('%s%s.log' % (self.Save_file_path, log_filename))
            stderr_streamHandler = logging.StreamHandler()

            # create formatter and add it to the handlers
            stderr_streamHandler.setFormatter(stderr_formatter)
            stderr_fileHandler.setFormatter(stderr_formatter)

            stderr_logger.addHandler(stderr_fileHandler)
            stderr_logger.addHandler(stderr_streamHandler)

            stderr_logger.setLevel(logging.ERROR)
            stderr_logger.error(e, exc_info)
            stderr_logger.error(message, exc_info)
        except EnvironmentError as err:
            stderr_logger = logging.getLogger()
            stderr_formatter = logging.Formatter('[%(process)d %(processName)s %(thread)d %(threadName)s|%(filename)s:%(lineno)s] %(asctime)s > %(message)s')
            stderr_fileHandler = logging.FileHandler('%s%s.log' % (self.Save_file_path, log_filename))
            stderr_streamHandler = logging.StreamHandler()

            # create formatter and add it to the handlers
            stderr_streamHandler.setFormatter(stderr_formatter)
            stderr_fileHandler.setFormatter(stderr_formatter)

            stderr_logger.addHandler(stderr_fileHandler)
            stderr_logger.addHandler(stderr_streamHandler)

            stderr_logger.setLevel(logging.ERROR)
            stderr_logger.error(err, exc_info)
            stderr_logger.error(message, exc_info)

class Inspector(object):
    def __init__(self, group, type="GS"):

        if type == "GS":
            self.groupName = "GS_" + group
            self.basedir = '/EQL7/pipeline/'
            self.bam_dir = '%s%s/' % (self.basedir, self.groupName)
            self.GS_test = "%sGS_test/" %(self.basedir)
            self.fastq_link = '/EQL2/%s/WXS/fastq/link/' %(self.groupName)
            self.BTscan = '/EQL2/btscan_pool/'
            self.GBMscan = '/EQL2/gbmscan_pool/'

        elif type == "RSq":
            self.groupName = "MACROGEN_" + group + "_rsq2"
            #self.groupName = "SGI" + group + "_rsq2"
            self.basedir = '/EQL8/pipeline/'
            self.tdf_dir = '%s%sexpr/' % (self.basedir, self.groupName)
            self.bam_dir = '%s%smut/' % (self.basedir, self.groupName)
            self.skip_dir = '%s%sskip/' % (self.basedir, self.groupName)
            self.ei_dir = '%s%seiJunc/' % (self.basedir, self.groupName)
            self.fusion_dir = '%s%sfusion/' % (self.basedir, self.groupName)
            self.fastq_link = '/EQL2/SGI_%s/RNASeq/fastq/link/' %(group)

        self.Pair = '_pair_filter_vep.dat'
        self.Single = '_single_filter_vep.dat'
        self.type = type

    def file_existance(self, fileLists, Dir):

        List_Done = []
        List_Not_Done = []

        for aa in range(len(fileLists)):
            #ReferenceDirs = glob(Dir+"*"+fileLists[aa]+"*" + self.Pair)
            ReferenceDirs = glob(Dir + "*" + fileLists[aa] + "*")
            if ReferenceDirs == []:
                if not fileLists[aa].split('_')[-1] == "B":
                    List_Not_Done.append(fileLists[aa])
            else:
                print(ReferenceDirs)
                List_Done.append(fileLists[aa])

        print('\n\n\n\n\n', List_Done, " exists!!!")
        print( '\n\n\n\n\n', List_Not_Done, " not exists")

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


    def Sample_naming_inspection(self, matchingList, isTAM=False):

        standards = []
        needToRename = []

        if self.type == "GS":
            standard1 = re.compile("^S([0-9]{2})([0-9]{5,})(_T)_(GS)+")
            standard2 = re.compile("^NS([0-9]{2})_([0-9]{3,4})(|_T0[1-9])(|_P[0-9]{1,2}|_S[0-9]|_SP[0-9]{1,2})_(GS)+")
            standard3 = re.compile("^NS([0-9]{2})_([0-9]{3,4})(|_T0[1-9])(_BR[1-9]|_BR[1-9]S[0-9])(|_P[0-9]{1,2}|_SP[0-9]{1,2})_(GS)+")
            standard4 = re.compile("^(|HGF_)(IRCR|SNU|CHA|NCC|AMC)_([A-Z]{2,3})([0-9]{2})_([0-9]{3,4})(_B|_T|_T0[1-9])_(GS)+")
            standard5 = re.compile("^(|HGF_)(IRCR|SNU|CHA|NCC|AMC)_([A-Z]{2,3})([0-9]{2})_([0-9]{3,4})(|_T|_T0[1-9])(|_P[0-9]{1,2}|_S[0-9]|_SP[0-9]{1,2})_(GS)+")
            standard6 = re.compile("^(|HGF_)(IRCR|SNU|CHA|NCC|AMC)_([A-Z]{2,3})([0-9]{2})_([0-9]{3,4})(|_T|_T0[1-9])(_BR[1-9]|_BR[1-9]S[0-9])(|_P[0-9]{1,2}|_SP[0-9]{1,2})_(GS)+")
            standards = [standard1, standard2, standard3, standard4, standard5, standard6]
        elif self.type == "RSq":
            standard1 = re.compile("^IRCR_([A-Z]{2,3})([0-9]{2})_([0-9]{3,4})(|_T0[1-9]|_T0[1-9]_P0)_(RSq)+")
            standard2 = re.compile("^IRCR_([A-Z]{2,3})([0-9]{2})_([0-9]{3,4})(_C0[0-9]|_T0[1-9]_C0[0-9])_(RSq)+") # Total cells, Tumor Cells
            standard3 = re.compile("^IRCR_([A-Z]{2,3})([0-9]{2})_([0-9]{3,4})(_TAM|_T0[1-9]_TAM)_(RSq)+") # TAM

            if isTAM == True:
                standards = [standard2, standard3]
            else:
                standards = [standard1]
        elif self.type == "WXS":
            standard1 = re.compile("^S([0-9]{2})([0-9]{5,})_T_(SS|TS)+")
            standard2 = re.compile("^NS([0-9]{2})_([0-9]{3,4})(|_T0[1-9])(|_P[0-9]{1,2}|_S[0-9]|_SP[0-9]{1,2})_(SS|TS)+")
            standard3 = re.compile("^NS([0-9]{2})_([0-9]{3,4})(|_T0[1-9])(_BR[1-9]|_BR[1-9]S[0-9])(|_P[0-9]{1,2}|_SP[0-9]{1,2})_(SS|TS)+")
            standard4 = re.compile("^(|HGF_)(IRCR|SNU|CHA|NCC|AMC)_([A-Z]{2,3})([0-9]{2})_([0-9]{3,4})(_B|_T|_T0[1-9])_(SS|TS)+")
            standard5 = re.compile("^(|HGF_)(IRCR|SNU|CHA|NCC|AMC)_([A-Z]{2,3})([0-9]{2})_([0-9]{3,4})(|_T0[1-9])(|_P[0-9]{1,2}|_S[0-9]|_SP[0-9]{1,2})_(SS|TS)+")
            standard6 = re.compile("^(|HGF_)(IRCR|SNU|CHA|NCC|AMC)_([A-Z]{2,3})([0-9]{2})_([0-9]{3,4})(|_T0[1-9])(_BR[1-9]|_BR[1-9]S[0-9])(|_P[0-9]{1,2}|_SP[0-9]{1,2})_(SS|TS)+")
            standards = [standard1, standard2, standard3, standard4, standard5, standard6]


        for ref in standards:
            needToRename = list(filter(lambda x: ref.match(x) == None, matchingList)) + needToRename

        return needToRename


    def Renaming_Files(self, isTAM=False):
        #sampleNames = glob("%s/*" % (self.bam_dir))
        #sampleNames = glob("%s/*" % (self.BTscan))
        #sampleNames = ["%sIRCR_S16_05210_GS" % (self.GS_test)]
        #sampleNames = glob("/EQL8/pipeline/TERAGEN_20160711_rsq2mut/*") + glob("/EQL8/pipeline/TERAGEN_20160801_rsq2mut/*") + glob("/EQL8/pipeline/TERAGEN_20160826_rsq2mut/*") + glob("/EQL8/pipeline/TERAGEN_20160906_rsq2mut/*") + glob("/EQL8/pipeline/TERAGEN_20160907_rsq2mut/*") + glob("/EQL8/pipeline/TERAGEN_20160912_rsq2mut/*")
        sampleNames = glob("/EQL8/pipeline/TERAGEN_20160801_rsq2expr/*") + \
                      glob("/EQL8/pipeline/TERAGEN_20160826_rsq2expr/*1076*") + glob("/EQL8/pipeline/TERAGEN_20160906_rsq2expr/*") + \
                      glob("/EQL8/pipeline/TERAGEN_20160907_rsq2expr/*") + glob("/EQL8/pipeline/TERAGEN_20160912_rsq2expr/*")


        #sampleNames = ["%sIRCR_BT16_1290_T02_P0_GS" %(self.GS_test)]
        unmatchedList = []

        ## At first, upper dir names are changed
        for sampleGroup in sampleNames:
            if sampleGroup.split(".")[-1] == "html":
                ######## Html #####################
                FirstFilter = sampleGroup.split("/")  # e.g : ['', 'EQL8', 'pipeline', 'SGI20170718', 'IRCR_BT16_1021_02_RSq.html']
                SecondFilter = FirstFilter[-1]  # e.g : ['IRCR_BT16_1021_02_RSq.html']
                HTMLFilter = self.Sample_naming_inspection([SecondFilter], True)  # e.g :
                HTMLFilter = list(set(HTMLFilter))
                if isTAM==True:
                    CorrectedHTML = self.correcting_naming_mistakes_TAM(HTMLFilter)
                else:
                    CorrectedHTML = self.correcting_naming_mistakes(HTMLFilter)

                ####### Link ######################
                if len(HTMLFilter) == 1:
                    OLD = "/".join(FirstFilter[:-1] + HTMLFilter)
                    NEW = "/".join(FirstFilter[:-1] + CorrectedHTML)

                    DIRFilter = HTMLFilter[0][:-5]
                    print(DIRFilter)
                    CorrectedName = CorrectedHTML[0][:-5]
                    #os.path.join(os.getcwd(), self.fastq_link)
                    OldLink1 = '%s%s.1.fq.gz' % (self.fastq_link, DIRFilter)
                    OldLink2 = '%s%s.2.fq.gz' % (self.fastq_link, DIRFilter)
                    NewLink1 = '%s%s.1.fq.gz' % (self.fastq_link, CorrectedName)
                    NewLink2 = '%s%s.2.fq.gz' % (self.fastq_link, CorrectedName)
                    print("old %s vs new %s" % (OldLink1, NewLink1))

                    self.Rename(OLD, NEW)
                    self.Rename(OldLink1, NewLink1)
                    self.Rename(OldLink2, NewLink2)
            else:
                continue

        ## Then, sub dir names and file names are changed
        for sampleGroup in sampleNames:
            matchingFiles = glob("%s/*" % (sampleGroup))  # e.g : ['/EQL8/pipeline/SGI20170718/IRCR_BT16_1021_02_RSq/IRCR_BT16_1021_02_RSq_splice.bam', '/EQL8/pipeline/SGI20170718/IRCR_BT16_1141_RSq/IRCR_BT16_1141_RSq_splice.dedup.bai']
            FirstFilter = list(map(methodcaller("split", "/"), matchingFiles)) # e.g : [['', 'EQL8', 'pipeline', 'SGI20170718', 'IRCR_BT16_1021_02_RSq', 'IRCR_BT16_1021_02_RSq_splice.bam'], ['', 'EQL8', 'pipeline', 'SGI20170718', 'IRCR_BT16_1141_RSq', 'IRCR_BT16_1141_RSq_splice.dedup.bai']]
            SecondFilter = list(map(lambda x: x[-1], FirstFilter)) # e.g : ['IRCR_BT16_1021_02_RSq_splice.bam', 'IRCR_BT16_1141_RSq_splice.dedup.bai']
            ThirdFilter = self.Sample_naming_inspection(SecondFilter, True) # e.g : ['IRCR_BT16_1021_02_RSq_splice.bam']
            ThirdFilter = list(set(ThirdFilter))
            if isTAM == True:
                CorretedFileNameList = self.correcting_naming_mistakes_TAM(ThirdFilter)
            else:
                CorretedFileNameList = self.correcting_naming_mistakes(ThirdFilter)
            OLDFILES = list(map(lambda x: sampleGroup + '/' + x, ThirdFilter))
            NEWFILES = list(map(lambda x: sampleGroup + '/' + x, CorretedFileNameList))

            if len(OLDFILES) == len(NEWFILES):
                for i in range(len(OLDFILES)):
                    self.Rename(OLDFILES[i], NEWFILES[i])

            FourthFilter = sampleGroup.split("/")
            FifthFilter = self.Sample_naming_inspection(FourthFilter[-1:], True)
            if isTAM == True:
                NEWDIRList = self.correcting_naming_mistakes_TAM(FifthFilter)
            else:
                NEWDIRList = self.correcting_naming_mistakes(FifthFilter)
            FourthFilter[-1] = NEWDIRList[-1] 
            self.Rename(sampleGroup, '/'.join(FourthFilter))



    def correcting_naming_mistakes_TAM(self, inputList):
        outputList = []
        mistake1 = re.compile("^IRCR_([A-Z]{2,3})([0-9]{2})_([0-9]{3,4})(|_T0[1-9]|_T0[1-9]_P0)_(RSq)+")

        for NAME in inputList:
            if mistake1.match(NAME) != None:
                newNAME = NAME.replace("_RSq", "_C00_RSq")
            else:
                newNAME = NAME
            outputList.append(newNAME)
        return outputList

    def correcting_naming_mistakes(self, inputList):
        outputList = []
        mistake1 = re.compile("^IRCR_([A-Z]{2,3})([0-9]{2})_([0-9]{3,4})(_0[1-9]|_0[1-9]_P0)_([A-Z]{2,3})+")
        mistake2 = re.compile("^S([0-9]{2})([0-9]{5,})_(GS)+")
        mistake3 = re.compile("^(|IRCR_)S([0-9]{2})_([0-9]{5,})_(GS)+")

        for NAME in inputList:
            if mistake1.match(NAME) != None:
                newNAME = NAME.replace("_0", "_T0")
            elif mistake2.match(NAME) != None:
                newNAME = NAME.replace("_GS", "_T_GS")
            elif mistake3.match(NAME) != None:
                if "IRCR_S" in NAME:
                    newNAME = NAME.replace("IRCR_S", "S")

                Delimiter = "_GS"
                newNAMEList = newNAME.split(Delimiter)
                newNAMEList[0] = ''.join(newNAMEList[0].split("_"))
                newNAME = Delimiter.join(newNAMEList)
                if mistake2.match(newNAME) != None:
                    newNAME = newNAME.replace("_GS", "_T_GS")
            else:
                newNAME = NAME
            outputList.append(newNAME)
        return outputList

    def ReLink(self,toUnlink, newLink, target=''):
        STDOUT = sp.check_output("ls -l %s" %(toUnlink), shell=True)
        if target == '':
            TargetPoint = str(STDOUT).split("->")[-1][1:-3]
        else:
            TargetPoint = target
            
        try:
            sp.check_call("unlink %s" % (toUnlink), shell=True)

        except FileNotFoundError as fnf:
            print("The symbolic link is already dead!!! %s" %(fnf))
            sp.check_call("rm -rf %s" % (toUnlink), shell=True)
            pass
        os.symlink(TargetPoint, newLink)


    def Rename(self, A, B):
        os.system('mv %s %s' % (A, B))
        print("%s >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> \n %s " %(A, B))

    def Manual_Renaming(self, target_Dir):
## self.GS_test >> self.bam_dir
        olddir = "%sIRCR_GBM14_682_BR1_P5_RSq" %(target_Dir)
        newdir = "%sIRCR_GBM15_682_BR1_P5_RSq" %(target_Dir)
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

    def read_Logs(self, oldExp, newExp):
        LogFileList = glob("/EQL8/pipeline/TERAGEN_20160711_rsq2expr/*C00_RSq/*log") + glob("/EQL8/pipeline/TERAGEN_20160801_rsq2expr/*C00_RSq/*log") + \
                      glob("/EQL8/pipeline/TERAGEN_20160826_rsq2expr/*1076*C00_RSq/*log") + glob("/EQL8/pipeline/TERAGEN_20160906_rsq2expr/*C00_RSq/*log") + \
                      glob("/EQL8/pipeline/TERAGEN_20160907_rsq2expr/*C00_RSq/*log") + glob("/EQL8/pipeline/TERAGEN_20160912_rsq2expr/*C00_RSq/*log")

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

###########################################################################################################
    def read_DAT(self, isPAIR=True):
        VEP = GS_variant_effect_prediction()
        VEP.initialize_setting(self.groupName)
        VEP.initialize_Single_Pair(isPair=isPAIR)
        LIST = VEP.Select_TODO_list()
        INDEX = ['Sample_Name', 'CHR', 'POSITION', 'REF', 'ALT', 'n_nRef', 'n_nAlt', 't_nRef', 't_nAlt', 'Annot_type', 'canonical', 'SC']
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
                            Sample = [SName, CHR, POS, REF, ALT, n_nRef, n_nAlt, t_nRef, t_nAlt, Annot_type, canonical, SC]
                            Sample_info = pd.Series(data=Sample, index=INDEX)
                            for element in Sample:
                                if element == '':
                                    print(line, Sample_info)
                                    instant_msg.append(file)
                                    with open("/home/jsgene/JK1/NGS/exec_history/%s_missing_in_DAT.txt" % ID, 'a') as wf:
                                        wf.write(line)
                                else:
                                    print(dirs, 'is OK')
                                    continue
                except IOError as err:
                    logging.basicConfig(filename="/home/jsgene/JK1/NGS/exec_history/%s_no_DAT_file.log" % ID, level=logging.DEBUG)
                    logging.debug(err)
                    continue
            print("something wrong with %s" %instant_msg)

#    def sequence_Quality(self):

class File_manager(object):
    import os, sys
    import subprocess
    def from_EQL_to_Storage(self, type):
        departureHost = "jsgene@119.86.100.106"
        departure = "/EQL8/pipeline/SGI20150709_rsq2expr/"
        destinationHost = "smcbi@119.4.212.148"
        destination = "/media/backup"
        if type == "GS":
            destination = destination + "/WXS_pipeline/"
        elif type == "RSq":
            destination = destination + "/RNASeq_pipeline/"

        cmd1 = 'sshpass -p smcbi '
        cmd2 = 'rsync -avXz --progress %s %s:%s' % (departure, destinationHost, destination)

        #cmd = cmd1 + cmd2
        cmd3 = "SSHPASS=smcbi"
        cmd4 = "sshpass -e sftp -oBatchMode=no -b - %s <<EOF" %(destinationHost)
        cmd5 = "put %s" %(departure)
        cmd6 = "bye"
        cmd7 = "EOF"
        print(cmd)
        os.system(cmd)


#################################################################################################
    def move_whole_dir(self, FROM, TO):
        sp.check_output(["cp" "-r", FROM, TO], stderr=sp.STDOUT, shell=True)

################################################################################################
    def Xeno_Link(self):
        mouse_fastqs = glob('/EQL2/*/WXS/fastq/*mouse*')
        link_exists = glob('/EQL2/*/WXS/fastq/link/*mm10_GS*')
        fastqs_hasLinks = {}
        for link in link_exists:
            destination = os.path.realpath(link)
            fastqs_hasLinks[destination] = link

        with open("/data1/home/jsgene/JSONS/Xeno_links.json", 'w') as f:
            json.dump(fastqs_hasLinks, f)

        Diction = pd.DataFrame(list(fastqs_hasLinks.items()), columns=['fastq', 'sym_link'])
        Diction.to_csv("/data1/home/jsgene/JSONS/Xeno_links.csv", sep=",")

        needTolink = []

        for fastq in mouse_fastqs:
            if fastq in fastqs_hasLinks.keys():
                continue
            else:
                if "mouse_1" in fastq:
                    originLink = fastq.replace("mouse_1", "human_1")
                elif "mouse_2" in fastq:
                    originLink = fastq.replace("mouse_2", "human_2")
                humanLink = originLink.split("/")[-1]
                filelist = originLink.split("WXS/fastq/")[0]
                needTolink.append((filelist, humanLink, fastq))

        if len(needTolink) > 0:
            for (superiorDir, humanLink, fastq) in needTolink:
                with open(superiorDir + "filelist.txt", 'r') as f:
                    for line in f:
                        print(line)
                        (ircr_name, DNA, sgilink) = line.split("\t")
                        if humanLink[:-12] in sgilink:
                            if "human_1" in humanLink:
                                symLink = superiorDir + "WXS/fastq/link/" + ircr_name + "_mm10_GS.1.fq.gz"
                            elif "human_2" in humanLink:
                                symLink = superiorDir + "WXS/fastq/link/" + ircr_name + "_mm10_GS.2.fq.gz"
                            try:
                                os.symlink(fastq, symLink)
                            except OSError as pe:
                                print (pe)
                                os.chmod(superiorDir, 0o777)
                                os.symlink(fastq, symLink)
                            print("%s <----- %s" % (fastq, symLink))

                        else:
                            print("There is no %s in %sfilelist.txt" % (humanLink[:-12], superiorDir))
        else:
            print("There is nothing to link")

#########################################################################################################################
    def Linking_multiple_dirs(self, actual_location=[]):
        for A_dir in actual_location:
            print(A_dir, " >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> ")
            loc_to_link = A_dir.replace("EQL11", "EQL7")
            print(loc_to_link)
            os.symlink(A_dir, loc_to_link)

    #####################################################################################################################

    def Link_modifier(self, prefix="SGI", query="C00", date = "20160826"):
        Origin = {"Fastq" : glob("/EQL2/%s_%s/RNASeq/fastq/link/*%s_RSq/*fq.gz" %(prefix, date, query)), "GSNAP": glob("/EQL8/pipeline/%s%s_rsq2mut/*%s_RSq/*splice.gsnap.gz" % (prefix, date, query))}
        ToUnlink = {"Fastq": glob("/EQL8/pipeline/%s_%s_rsq2mut/*%s_RSq/*fq.gz" %(prefix, date, query)),
                    "GSNAP": glob("/EQL8/pipeline/%s_%s_rsq2skip/*%s_RSq/*splice.gsnap.gz" %(prefix, date, query)) +
                             glob("/EQL8/pipeline/%s_%s_rsq2eiJunc/*%s_RSq/*splice.gsnap.gz" %(prefix, date, query)) +
                             glob("/EQL8/pipeline/%s_%s_rsq2fusion/*%s_RSq/*splice.gsnap.gz" % (prefix, date, query))}

        for unlink in ToUnlink["Fastq"]:
            print(unlink)
            target_matching = unlink.split("/")[-1]
            for origin in Origin:
                target = origin.split("/")[-1]
                if target == target_matching:
                    ins = Inspector("20170830", "RSq")
                    ins.ReLink(unlink, unlink, origin)

    def managing_JSONs(self):
        # samples = ["IRCR_BT16_1061_B", "IRCR_BT16_1061_T", "IRCR_BT16_943_B", "IRCR_BT16_943_T", "IRCR_BT15_842_B", "IRCR_BT15_842_T", "IRCR_BT16_1165_B", "IRCR_BT16_1165_T", "IRCR_GBM12_140_BR1", "IRCR_BT16_966_B", "IRCR_BT16_966_T02", "IRCR_BT16_966_T04", "IRCR_GBM15_717_B", "IRCR_GBM15_717_T", "IRCR_BT16_936_B", "IRCR_BT16_936_T01", "IRCR_BT16_936_T02", "IRCR_BT15_884_B", "IRCR_BT15_884_T01", "IRCR_BT15_884_T05", "IRCR_BT16_1047_B", "IRCR_BT16_1047_T", "IRCR_BT15_890_B", "IRCR_BT15_890_T", "IRCR_BT16_937_B", "IRCR_BT16_937_T", "IRCR_BT16_1045_B", "IRCR_BT16_1045_T01", "IRCR_BT16_1045_T02", "IRCR_BT16_1045_T03", "IRCR_BT16_1045_T05", "IRCR_BT16_1072_B", "IRCR_BT16_1072_T01", "IRCR_BT16_1072_T02", "IRCR_BT16_1081_B", "IRCR_BT16_1081_T02", "IRCR_BT16_1105_B", "IRCR_BT16_1105_T01", "IRCR_BT16_1105_T02", "IRCR_BT16_948_B", "IRCR_BT16_948_T01", "IRCR_BT16_948_T02", "IRCR_BT16_1112_B", "IRCR_BT16_1112_T", "IRCR_BT16_945_B", "IRCR_BT16_945_T", "IRCR_BT16_951_B", "IRCR_BT16_951_T01", "IRCR_BT16_951_T03", "IRCR_BT16_951_T05", "IRCR_BT16_1059_B", "IRCR_BT16_1059_T01", "IRCR_BT16_1059_T02", "IRCR_BT16_1078_B", "IRCR_BT16_1078_T", "IRCR_BT16_1141_B", "IRCR_BT16_1141_T", "IRCR_BT16_1144_B", "IRCR_BT16_1144_T", "IRCR_BT16_1146_B", "IRCR_BT16_1146_T", "IRCR_BT16_1147_B", "IRCR_BT16_1147_T", "IRCR_BT16_1149_B", "IRCR_BT16_1149_T", "IRCR_BT16_1153_B", "IRCR_BT16_1153_T01", "IRCR_BT16_1153_T02", "IRCR_BT16_1153_T03", "IRCR_BT16_1160_B", "IRCR_BT16_1160_T01", "IRCR_BT16_1160_T02", "IRCR_BT16_1170_B", "IRCR_BT16_1170_T", "IRCR_BT16_1173_B", "IRCR_BT16_1173_T01", "IRCR_BT16_1173_T02", "IRCR_BT16_1174_B", "IRCR_BT16_1174_T", "IRCR_BT16_1179_B", "IRCR_BT16_1179_T", "IRCR_BT16_1185_B", "IRCR_BT16_1185_T", "IRCR_BT16_1197_B", "IRCR_BT16_1197_T", "IRCR_BT16_1199_B", "IRCR_BT16_1199_T", "IRCR_BT16_995_B", "IRCR_BT16_995_T", "IRCR_BT16_1140_B", "IRCR_BT16_1140_T", "IRCR_BT16_1143_B", "IRCR_BT16_1143_T", "IRCR_BT16_1217_B", "IRCR_BT16_1217_T", "IRCR_BT16_1205_B", "IRCR_BT16_1205_T01", "IRCR_BT16_1184_B", "IRCR_BT16_1184_T", "IRCR_BT16_1195_B", "IRCR_BT16_1195_T", "IRCR_BT16_1193_B", "IRCR_BT16_1193_T", "IRCR_BT16_1186_B", "IRCR_BT16_1186_T01", "IRCR_BT16_1186_T02", "IRCR_BT16_1208_B", "IRCR_BT16_1208_T", "IRCR_BT15_904_B", "IRCR_BT15_904_T01", "IRCR_BT15_904_T04", "IRCR_BT15_904_T05", "IRCR_BT16_1218_B", "IRCR_BT16_1218_T", "IRCR_BT16_1231_B", "IRCR_BT16_1231_T", "IRCR_BT16_1232_B", "IRCR_BT16_1232_T01", "IRCR_BT16_1232_T04", "IRCR_BT16_1236_B", "IRCR_BT16_1236_T", "IRCR_BT16_1237_B", "IRCR_BT16_1237_T", "IRCR_BT16_1240_B", "IRCR_BT16_1240_T", "IRCR_BT16_1245_B", "IRCR_BT16_1245_T02", "IRCR_BT16_1249_B", "IRCR_BT16_1249_T01", "IRCR_BT16_1249_T02", "IRCR_BT16_1250_B", "IRCR_BT16_1250_T", "IRCR_BT17_1258_B", "IRCR_BT17_1258_T02", "NCC_BT17_001_B", "NCC_BT17_001_T03", "NCC_BT17_001_T02", "NCC_BT17_001_T01", "IRCR_BT17_1269_B", "IRCR_BT17_1269_T01", "IRCR_BT17_1269_T02", "IRCR_BT17_1271_B", "IRCR_BT17_1271_T", "IRCR_BT17_1279_B", "IRCR_BT17_1279_T", "IRCR_BT17_1261_B", "IRCR_BT17_1261_T", "IRCR_BT17_1262_B", "IRCR_BT17_1262_T", "IRCR_BT17_1264_B", "IRCR_BT17_1264_T", "IRCR_BT17_1267_B", "IRCR_BT17_1267_T", "IRCR_BT17_1277_B", "IRCR_BT17_1277_T", "IRCR_BT17_1287_B", "IRCR_BT17_1287_T", "IRCR_BT17_1290_B", "IRCR_BT17_1290_T04", "IRCR_BT17_1290_T05", "IRCR_BT17_1290_T06", "IRCR_BT17_1291_B", "IRCR_BT17_1291_T02", "IRCR_BT17_1291_T03", "IRCR_BT17_1297_B", "IRCR_BT17_1297_T", "IRCR_BT17_1299_B", "IRCR_BT17_1299_T", "IRCR_BT17_1301_B", "IRCR_BT17_1301_T", "IRCR_BT17_1305_B", "IRCR_BT17_1305_T02", "IRCR_BT17_1305_T03", "IRCR_BT17_1305_T04", "CHA_BT17_002_B", "CHA_BT17_002_T", "SNU_BT17_003_B", "SNU_BT17_003_T", "SNU_BT17_008_B", "SNU_BT17_008_T", "SNU_BT17_009_B", "SNU_BT17_009_T", "NCC_BT17_011_B", "NCC_BT17_011_T02", "AMC_BT17_005_B", "AMC_BT17_005_T01", "AMC_BT17_006_B", "AMC_BT17_006_T", "AMC_BT17_007_B", "AMC_BT17_007_T", "IRCR_BT17_1303_B", "IRCR_BT17_1303_T02", "IRCR_BT17_1303_T03", "IRCR_BT17_1306_B", "IRCR_BT17_1306_T", "IRCR_BT17_1307_B", "IRCR_BT17_1307_T", "IRCR_BT15_891_B", "IRCR_BT15_891_T01", "SNU_BT17_020_B", "SNU_BT17_020_T", "IRCR_BT17_1324_B", "IRCR_BT17_1324_T", "IRCR_BT17_1326_B", "IRCR_BT17_1326_T", "IRCR_BT17_1330_B", "IRCR_BT17_1330_T01", "IRCR_BT17_1330_T02", "IRCR_BT17_1330_T03", "IRCR_BT17_1330_T06", "IRCR_BT17_1332_B", "IRCR_BT17_1332_T", "LMJ_VIRUS_002_T", "IRCR_BT17_1317_B", "IRCR_BT17_1317_T01", "IRCR_BT17_1331_B", "IRCR_BT17_1331_T04", "IRCR_BT17_1339_B", "IRCR_BT17_1339_T01", "IRCR_BT17_1339_T02", "IRCR_BT17_1339_T03", "IRCR_BT17_1339_T04", "CHA_BT17_016_B", "CHA_BT17_016_T01", "CHA_BT17_016_T02", "CHA_BT17_016_T03", "CHA_BT17_016_T04", "SNU_BT17_013_B", "SNU_BT17_013_T", "SNU_BT17_014_B", "SNU_BT17_014_T", "SNU_BT17_015_B", "SNU_BT17_015_T", "SNU_BT17_017_B", "SNU_BT17_017_T", "IRCR_BT17_1272_B", "IRCR_BT17_1272_T", "IRCR_GBM14_400_B", "IRCR_GBM14_400_T", "IRCR_GBM15_766_B", "IRCR_GBM15_766_T", "IRCR_GBM14_608_BR1S2", "IRCR_BT17_1332_B", "IRCR_BT17_1332_T", "IRCR_GBM13_316_BR1S1", "IRCR_GBM15_698_BR1S2", "IRCR_BT15_849_T02_BR1S1", "IRCR_BT16_930_T02_BR1S1", "IRCR_GBM15_677_BR1S2", "IRCR_BT15_802_T03_BR1S1", "IRCR_GBM14_542_T01_BR1S1", "IRCR_BT15_891_T01_BR1S2", "IRCR_BT16_1134_T", "IRCR_BT16_1154_T", "IRCR_BT15_1134_B", "IRCR_BT16_1154_B", "IRCR_BT17_1290_T01", "IRCR_BT17_1342_B", "IRCR_BT17_1342_T", "IRCR_BT17_1343_B", "IRCR_BT17_1343_T", "IRCR_BT17_1346_B", "IRCR_BT17_1346_T01", "IRCR_BT17_1350_B", "IRCR_BT17_1350_T", "IRCR_BT17_1354_B", "IRCR_BT17_1354_T01", "IRCR_BT17_1354_T02", "IRCR_BT17_1357_B", "IRCR_BT17_1357_T", "SNU_BT17_022_B", "SNU_BT17_022_T01", "SNU_BT17_022_T02", "AMC_BT17_021_B", "AMC_BT17_021_T01", "AMC_BT17_021_T02", "IRCR_BT17_1303_T01", "IRCR_BT17_1353_B", "IRCR_BT17_1353_T01", "IRCR_BT17_1353_T02", "IRCR_BT17_1353_T03", "IRCR_BT16_1173_T06_P0", "IRCR_BT16_1245_T01_P0", "IRCR_BT17_1247_B", "IRCR_BT16_1247_T", "IRCR_BT16_1290_T02_P0", "IRCR_BT16_1177_T", "IRCR_BT17_1358_B", "IRCR_BT17_1358_T"]
        samples = ["S1318157_T", "S1605210_T", "IRCR_BT17_1319_B", "IRCR_BT15_870_B", "IRCR_BT15_870_T",
                   "IRCR_BT16_1023_B", "IRCR_BT16_1023_T", "IRCR_GBM13_300_B", "IRCR_GBM13_300_T", "IRCR_GBM14_509_B",
                   "IRCR_GBM14_509_T", "IRCR_GBM14_510_B", "IRCR_GBM14_510_T", "IRCR_BT16_935_B", "IRCR_BT16_935_T",
                   "IRCR_BT15_925_B", "IRCR_BT15_925_T", "IRCR_BT17_1323_B", "IRCR_BT17_1323_T01", "IRCR_BT17_1323_T02",
                   "IRCR_BT17_1323_T04", "IRCR_BT17_1325_B", "IRCR_BT17_1325_T01", "IRCR_BT17_1325_T03",
                   "IRCR_GBM10_038_B", "IRCR_GBM10_038_T", "IRCR_GBM13_231_B", "IRCR_GBM13_231_T", "IRCR_GBM13_292_B",
                   "IRCR_GBM13_292_T", "IRCR_GBM13_319_B", "IRCR_GBM13_319_T", "IRCR_GBM15_711_B", "IRCR_GBM15_711_T",
                   "IRCR_GBM14_458_B", "IRCR_GBM14_458_T", "IRCR_GBM14_487_B", "IRCR_GBM14_487_T", "IRCR_GBM14_494_B",
                   "IRCR_GBM14_494_T", "IRCR_GBM14_497_B", "IRCR_GBM14_497_T", "IRCR_GBM14_508_B", "IRCR_GBM14_508_T",
                   "IRCR_GBM14_511_B", "IRCR_GBM14_511_T", "IRCR_GBM14_517_B", "IRCR_GBM14_517_T02", "IRCR_GBM14_571_B",
                   "IRCR_GBM14_571_T", "IRCR_GBM14_614_B", "IRCR_GBM14_614_T", "IRCR_GBM14_616_B", "IRCR_GBM14_616_T",
                   "IRCR_GBM14_618_B", "IRCR_GBM14_618_T", "IRCR_GBM15_686_B", "IRCR_GBM15_686_T", "IRCR_GBM14_534_B",
                   "IRCR_GBM14_534_T", "IRCR_GBM14_589_B", "IRCR_GBM14_589_T", "IRCR_GBM14_567_B", "IRCR_GBM14_567_T",
                   "IRCR_GBM14_586_B", "IRCR_GBM14_586_T", "IRCR_GBM14_626_B", "IRCR_GBM14_626_T", "CHA_BT17_023_B",
                   "CHA_BT17_023_T03", "CHA_BT17_023_T04", "SNU_BT17_026_B", "SNU_BT17_026_T01", "SNU_BT17_027_B",
                   "SNU_BT17_027_T", "IRCR_BT17_1379_B", "IRCR_BT17_1379_T", "IRCR_BT17_1384_B", "IRCR_BT17_1384_T",
                   "IRCR_BT17_1351_B", "S16124067_T", "SNU_BT17_024_B", "SNU_BT17_024_T01", "SNU_BT17_025_B",
                   "SNU_BT17_025_T01", "IRCR_BT17_1366_B", "IRCR_BT17_1366_T", "IRCR_BT17_1367_B", "IRCR_BT17_1367_T",
                   "IRCR_BT17_1373_B", "IRCR_BT17_1373_T01", "IRCR_BT17_1373_T02", "IRCR_BT17_1373_T04",
                   "IRCR_BT17_1375_B", "IRCR_BT17_1375_T", "IRCR_BT17_1340_B", "IRCR_BT17_1340_T"]
        unknown_sample = []
        for i, name in enumerate(samples):

            temp_dict = {}
            path = "/data1/home/jsgene/JSONS/%s.json" % name
            try:
                with open(path, "r") as rf:
                    temp_dict = json.load(rf)
                    oldnum = temp_dict["No"]
                    temp_dict["No"] = i + 1
                    print("%s : number has just been changed from %d to %d" % (name, oldnum, i + 1))
                    temp_dict["Scan_Type"] = "BT_Scan"
                with open(path, "w") as wf:
                    json.dump(temp_dict, wf)
                    # os.chmod(path, 0o777)
            except FileNotFoundError as fnf:
                print(fnf)
                unknown_sample.append(path)




########################################################################################################################

if __name__=="__main__":

    ############################ Change Dir names and File names #######################################################
    ##ins = Inspector("20170516", "GS")
    #ins = Inspector("20170830", "RSq")
    #ins.Renaming_Files(isTAM=True)
    #ins.read_Logs("_RSq", "_C00_RSq")  ## Editing Log Files
    ###ins.Manual_Renaming(ins.bam_dir)

    ############################# Reset Symbolic Links and Create again ################################################
    # ins.ReLink('/EQL7/pipeline/GS_test/S16124067_T_GS/S16124067_T_GS.recal.bai', '/EQL7/pipeline/GS_test/S16124067_T_GS/S16124067_T_GS.recal.bai', '/EQL7/pipeline/GS_20170711/S16124067_T_GS/S16124067_T_GS.recal.bai')


    ########################### Find out duplicated dirs ###############################################################
    # ins.Duplicated_Dirs(["IRCR_BT16_1061_B", "IRCR_BT16_1061_T", "IRCR_BT16_943_B", "IRCR_BT16_943_T", "IRCR_BT15_842_B", "IRCR_BT15_842_T", "IRCR_BT16_1165_B", "IRCR_BT16_1165_T", "IRCR_GBM12_140_BR1", "IRCR_BT16_966_B", "IRCR_BT16_966_T02", "IRCR_BT16_966_T04", "IRCR_GBM15_717_B", "IRCR_GBM15_717_T", "IRCR_BT16_936_B", "IRCR_BT16_936_T01", "IRCR_BT16_936_T02", "IRCR_BT15_884_B", "IRCR_BT15_884_T01", "IRCR_BT15_884_T05", "IRCR_BT16_1047_B", "IRCR_BT16_1047_T", "IRCR_BT15_890_B", "IRCR_BT15_890_T", "IRCR_BT16_937_B", "IRCR_BT16_937_T", "IRCR_BT16_1045_B", "IRCR_BT16_1045_T01", "IRCR_BT16_1045_T02", "IRCR_BT16_1045_T03", "IRCR_BT16_1045_T05", "IRCR_BT16_1072_B", "IRCR_BT16_1072_T01", "IRCR_BT16_1072_T02", "IRCR_BT16_1081_B", "IRCR_BT16_1081_T02", "IRCR_BT16_1105_B", "IRCR_BT16_1105_T01", "IRCR_BT16_1105_T02", "IRCR_BT16_948_B", "IRCR_BT16_948_T01", "IRCR_BT16_948_T02", "IRCR_BT16_1112_B", "IRCR_BT16_1112_T", "IRCR_BT16_945_B", "IRCR_BT16_945_T", "IRCR_BT16_951_B", "IRCR_BT16_951_T01", "IRCR_BT16_951_T03", "IRCR_BT16_951_T05", "IRCR_BT16_1059_B", "IRCR_BT16_1059_T01", "IRCR_BT16_1059_T02", "IRCR_BT16_1078_B", "IRCR_BT16_1078_T", "IRCR_BT16_1141_B", "IRCR_BT16_1141_T", "IRCR_BT16_1144_B", "IRCR_BT16_1144_T", "IRCR_BT16_1146_B", "IRCR_BT16_1146_T", "IRCR_BT16_1147_B", "IRCR_BT16_1147_T", "IRCR_BT16_1149_B", "IRCR_BT16_1149_T", "IRCR_BT16_1153_B", "IRCR_BT16_1153_T01", "IRCR_BT16_1153_T02", "IRCR_BT16_1153_T03", "IRCR_BT16_1160_B", "IRCR_BT16_1160_T01", "IRCR_BT16_1160_T02", "IRCR_BT16_1170_B", "IRCR_BT16_1170_T", "IRCR_BT16_1173_B", "IRCR_BT16_1173_T01", "IRCR_BT16_1173_T02", "IRCR_BT16_1174_B", "IRCR_BT16_1174_T", "IRCR_BT16_1179_B", "IRCR_BT16_1179_T", "IRCR_BT16_1185_B", "IRCR_BT16_1185_T", "IRCR_BT16_1197_B", "IRCR_BT16_1197_T", "IRCR_BT16_1199_B", "IRCR_BT16_1199_T", "IRCR_BT16_995_B", "IRCR_BT16_995_T", "IRCR_BT16_1140_B", "IRCR_BT16_1140_T", "IRCR_BT16_1143_B", "IRCR_BT16_1143_T", "IRCR_BT16_1217_B", "IRCR_BT16_1217_T", "IRCR_BT16_1205_B", "IRCR_BT16_1205_T01", "IRCR_BT16_1184_B", "IRCR_BT16_1184_T", "IRCR_BT16_1195_B", "IRCR_BT16_1195_T", "IRCR_BT16_1193_B", "IRCR_BT16_1193_T", "IRCR_BT16_1186_B", "IRCR_BT16_1186_T01", "IRCR_BT16_1186_T02", "IRCR_BT16_1208_B", "IRCR_BT16_1208_T", "IRCR_BT15_904_B", "IRCR_BT15_904_T01", "IRCR_BT15_904_T04", "IRCR_BT15_904_T05", "IRCR_BT16_1218_B", "IRCR_BT16_1218_T", "IRCR_BT16_1231_B", "IRCR_BT16_1231_T", "IRCR_BT16_1232_B", "IRCR_BT16_1232_T01", "IRCR_BT16_1232_T04", "IRCR_BT16_1236_B", "IRCR_BT16_1236_T", "IRCR_BT16_1237_B", "IRCR_BT16_1237_T", "IRCR_BT16_1240_B", "IRCR_BT16_1240_T", "IRCR_BT16_1245_B", "IRCR_BT16_1245_T02", "IRCR_BT16_1249_B", "IRCR_BT16_1249_T01", "IRCR_BT16_1249_T02", "IRCR_BT16_1250_B", "IRCR_BT16_1250_T", "IRCR_BT17_1258_B", "IRCR_BT17_1258_T02", "NCC_BT17_001_B", "NCC_BT17_001_T03", "NCC_BT17_001_T02", "NCC_BT17_001_T01", "IRCR_BT17_1269_B", "IRCR_BT17_1269_T01", "IRCR_BT17_1269_T02", "IRCR_BT17_1271_B", "IRCR_BT17_1271_T", "IRCR_BT17_1279_B", "IRCR_BT17_1279_T", "IRCR_BT17_1261_B", "IRCR_BT17_1261_T", "IRCR_BT17_1262_B", "IRCR_BT17_1262_T", "IRCR_BT17_1264_B", "IRCR_BT17_1264_T", "IRCR_BT17_1267_B", "IRCR_BT17_1267_T", "IRCR_BT17_1277_B", "IRCR_BT17_1277_T", "IRCR_BT17_1287_B", "IRCR_BT17_1287_T", "IRCR_BT17_1290_B", "IRCR_BT17_1290_T04", "IRCR_BT17_1290_T05", "IRCR_BT17_1290_T06", "IRCR_BT17_1291_B", "IRCR_BT17_1291_T02", "IRCR_BT17_1291_T03", "IRCR_BT17_1297_B", "IRCR_BT17_1297_T", "IRCR_BT17_1299_B", "IRCR_BT17_1299_T", "IRCR_BT17_1301_B", "IRCR_BT17_1301_T", "IRCR_BT17_1305_B", "IRCR_BT17_1305_T02", "IRCR_BT17_1305_T03", "IRCR_BT17_1305_T04", "CHA_BT17_002_B", "CHA_BT17_002_T", "SNU_BT17_003_B", "SNU_BT17_003_T", "SNU_BT17_008_B", "SNU_BT17_008_T", "SNU_BT17_009_B", "SNU_BT17_009_T", "NCC_BT17_011_B", "NCC_BT17_011_T02", "AMC_BT17_005_B", "AMC_BT17_005_T01", "AMC_BT17_006_B", "AMC_BT17_006_T", "AMC_BT17_007_B", "AMC_BT17_007_T", "IRCR_BT17_1303_B", "IRCR_BT17_1303_T02", "IRCR_BT17_1303_T03", "IRCR_BT17_1306_B", "IRCR_BT17_1306_T", "IRCR_BT17_1307_B", "IRCR_BT17_1307_T", "IRCR_BT15_891_B", "IRCR_BT15_891_T01", "IRCR_S13_18157", "IRCR_S16_05210", "IRCR_BT17_1319_B", "SNU_BT17_020_B", "SNU_BT17_020_T", "IRCR_BT17_1324_B", "IRCR_BT17_1324_T", "IRCR_BT17_1326_B", "IRCR_BT17_1326_T", "IRCR_BT17_1330_B", "IRCR_BT17_1330_T01", "IRCR_BT17_1330_T02", "IRCR_BT17_1330_T03", "IRCR_BT17_1330_T06", "IRCR_BT17_1332_B", "IRCR_BT17_1332_T", "IRCR_BT15_870_B", "IRCR_BT15_870_T", "IRCR_BT16_1023_B", "IRCR_BT16_1023_T", "IRCR_GBM13_300_B", "IRCR_GBM13_300_T", "LMJ_VIRUS_002_T", "IRCR_BT17_1317_B", "IRCR_BT17_1317_T01", "IRCR_BT17_1331_B", "IRCR_BT17_1331_T04", "IRCR_BT17_1339_B", "IRCR_BT17_1339_T01", "IRCR_BT17_1339_T02", "IRCR_BT17_1339_T03", "IRCR_BT17_1339_T04", "CHA_BT17_016_B", "CHA_BT17_016_T01", "CHA_BT17_016_T02", "CHA_BT17_016_T03", "CHA_BT17_016_T04", "SNU_BT17_013_B", "SNU_BT17_013_T", "SNU_BT17_014_B", "SNU_BT17_014_T", "SNU_BT17_015_B", "SNU_BT17_015_T", "SNU_BT17_017_B", "SNU_BT17_017_T", "IRCR_BT17_1272_B", "IRCR_BT17_1272_T", "IRCR_GBM14_400_B", "IRCR_GBM14_400_T", "IRCR_GBM15_766_B", "IRCR_GBM15_766_T", "IRCR_GBM14_608_BR1S2", "IRCR_BT16_935_B", "IRCR_BT16_935_T", "IRCR_BT15_925_B", "IRCR_BT15_925_T", "IRCR_BT17_1323_B", "IRCR_BT17_1323_T01", "IRCR_BT17_1323_T02", "IRCR_BT17_1323_T04", "IRCR_BT17_1325_B", "IRCR_BT17_1325_T01", "IRCR_BT17_1325_T03", "IRCR_BT17_1332_B", "IRCR_BT17_1332_T", "IRCR_GBM10_038_B", "IRCR_GBM10_038_T", "IRCR_GBM13_231_B", "IRCR_GBM13_231_T", "IRCR_GBM13_292_B", "IRCR_GBM13_292_T", "IRCR_GBM13_319_B", "IRCR_GBM13_319_T", "IRCR_GBM15_711_B", "IRCR_GBM15_711_T", "IRCR_GBM14_458_B", "IRCR_GBM14_458_T", "IRCR_GBM14_487_B", "IRCR_GBM14_487_T", "IRCR_GBM14_494_B", "IRCR_GBM14_494_T", "IRCR_GBM14_497_B", "IRCR_GBM14_497_T", "IRCR_GBM14_508_B", "IRCR_GBM14_508_T", "IRCR_GBM13_316_BR1S1", "IRCR_GBM15_698_BR1S2", "IRCR_BT15_849_T02_BR1S1", "IRCR_BT16_930_T02_BR1S1", "IRCR_GBM15_677_BR1S2", "IRCR_BT15_802_T03_BR1S1", "IRCR_GBM14_542_T01_BR1S1", "IRCR_BT15_891_T01_BR1S2", "IRCR_GBM14_589_B", "IRCR_GBM14_589_T", "IRCR_GBM14_567_B", "IRCR_GBM14_567_T", "IRCR_GBM14_586_B", "IRCR_GBM14_586_T", "IRCR_GBM14_626_B", "IRCR_GBM14_626_T", "IRCR_BT16_1134_T", "IRCR_BT16_1154_T", "IRCR_BT15_1134_B", "IRCR_BT16_1154_B", "IRCR_BT17_1290_T01", "IRCR_BT17_1342_B", "IRCR_BT17_1342_T", "IRCR_BT17_1343_B", "IRCR_BT17_1343_T", "IRCR_BT17_1346_B", "IRCR_BT17_1346_T01", "IRCR_BT17_1350_B", "IRCR_BT17_1350_T", "IRCR_BT17_1354_B", "IRCR_BT17_1354_T01", "IRCR_BT17_1354_T02", "IRCR_BT17_1357_B", "IRCR_BT17_1357_T", "SNU_BT17_022_B", "SNU_BT17_022_T01", "SNU_BT17_022_T02", "AMC_BT17_021_B", "AMC_BT17_021_T01", "AMC_BT17_021_T02", "IRCR_BT17_1303_T01", "IRCR_BT17_1353_B", "IRCR_BT17_1353_T01", "IRCR_BT17_1353_T02", "IRCR_BT17_1353_T03", "CHA_BT17_023_B", "CHA_BT17_023_T03", "CHA_BT17_023_T04", "SNU_BT17_026_B", "SNU_BT17_026_T01", "SNU_BT17_027_B", "SNU_BT17_027_T", "IRCR_BT17_1379_B", "IRCR_BT17_1379_T", "IRCR_BT17_1384_B", "IRCR_BT17_1384_T", "IRCR_BT17_1351_B", "S16124067_T", "SNU_BT17_024_B", "SNU_BT17_024_T01", "SNU_BT17_025_B", "SNU_BT17_025_T01", "IRCR_BT17_1366_B", "IRCR_BT17_1366_T", "IRCR_BT17_1367_B", "IRCR_BT17_1367_T", "IRCR_BT17_1373_B", "IRCR_BT17_1373_T01", "IRCR_BT17_1373_T02", "IRCR_BT17_1373_T04", "IRCR_BT17_1375_B", "IRCR_BT17_1375_T", "IRCR_BT17_1340_B", "IRCR_BT17_1340_T", "IRCR_BT16_1173_T06_P0", "IRCR_BT16_1245_T01_P0", "IRCR_BT17_1247_B", "IRCR_BT16_1247_T", "IRCR_BT16_1290_T02_P0", "IRCR_BT16_1177_T"])

    ########################### Inspect NGS output files whether their forms and dimensions interrupted or not #########
    # ins.read_DAT()

    ###########################
    mgr = File_manager(glob("/EQL11/pipeline/*"))
    ########################## Move Files to IRCR storage server #######################################################
    # mgr.from_EQL_to_Storage("RSq")

    ############################ Create Multiple Symbolic Links ########################################################
    mgr.Linking_multiple_dirs()

    ############################# Reset Symbolic Links and Create again ################################################
    #mgr.Link_modifier("TERAGEN", "C00", "20160826")


