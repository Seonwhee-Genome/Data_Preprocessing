import os, re, sys, stat
import json, pickle
import subprocess as sp
from glob import glob
from operator import methodcaller

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
        self.Link_Memo_isoforms = "/home/jsgene/JK1/NGS/20171124_link1.pkl"
        self.Link_Memo_gsnap = "/home/jsgene/JK1/NGS/20171124_link2.pkl"
        self.Link_Memo_fastq = "/home/jsgene/JK1/NGS/20171124_link3.pkl"

    def Rename_Rootdirectories_EQL2(self, sampleNum, replaceFrom, replaceTo):
        TargetList = glob(self.fastq_link+"*%s*"%sampleNum)
        for target in TargetList:
            self.Rename(target, target.replace(replaceFrom, replaceTo))

    def Rename_Subdirectories_EQL8(self, sampleNum, replaceFrom, replaceTo, Target, STEP):
        if Target == "Files":
            option = "*%s*/*" % sampleNum
        elif Target == "Dir":
            option = "*%s*" % sampleNum

        if STEP == 1:
            subList = glob(self.ei_dir+option) + glob(self.fusion_dir+option) + glob(self.skip_dir+option)
        elif STEP == 2:
            subList = glob(self.bam_dir+option)
        elif STEP == 3:
            subList = glob(self.tdf_dir + option)

        FileShouldBeRenamed = []
        for paths in subList:
            Splited = paths.split("/")
            changing_word = Splited[-1]
            roots = "/".join(Splited[:-1]) + "/"
            changed_word = changing_word.replace(replaceFrom, replaceTo)
            """path에서 /를 기준으로 가장 끝단부터 이름을 바꿔줘야 하기 때문에   (예:  /폴더명/파일명 --> /폴더명/바꾼파일명 --> /바꾼폴더명/바꾼파일명)
            이름 바꿀 부분과 바꾸지 않고 보존할 부분을 구분하는 것이다"""

            FileShouldBeRenamed.append((roots+changing_word, roots+changed_word))

        for a_tuple in FileShouldBeRenamed:
            self.Rename(a_tuple[0], a_tuple[1])

    def making_Unlink(self, query, STEP):
        if STEP == 1:
            OldLinks = glob("%s*%s*/*splice.gsnap.gz" % (self.skip_dir, query)) + glob("%s*%s*/*splice.gsnap.gz" % (self.ei_dir, query)) + glob("%s*%s*/*splice.gsnap.gz" % (self.fusion_dir, query))
            Link_Memo = self.Link_Memo_isoforms
        elif STEP == 2:
            OldLinks = glob("%s*%s*/*fq.gz" % (self.bam_dir, query))
            Link_Memo = self.Link_Memo_gsnap
        elif STEP == 3:
            OldLinks = glob("%s*%s*/*fq.gz" % (self.tdf_dir, query))
            Link_Memo = self.Link_Memo_fastq
        with open(Link_Memo, "wb") as f:
            pickle.dump(OldLinks, f)
        for LinkToBeRemoved in OldLinks:
            try:
                print("unlink ", LinkToBeRemoved)
                sp.check_call("unlink %s" % (LinkToBeRemoved), shell=True)

            except FileNotFoundError as fnf:
                print("The symbolic link is already dead!!! %s" % (fnf))
                sp.check_call("rm -rf %s" % (LinkToBeRemoved), shell=True)
                pass

    def construct_Link(self, query, STEP):
        Origin = {"Fastq": glob("%s/*%s*fq.gz" % (self.fastq_link, query)),
                  "GSNAP": glob("%s*%s*_RSq/*splice.gsnap.gz" % (self.bam_dir, query))}
        if STEP == 4:
            key = "GSNAP"
            Link_Memo = self.Link_Memo_isoforms
            Link_source = Origin[key]
        elif STEP == 5:
            key = "Fastq"
            Link_Memo = self.Link_Memo_fastq
            Link_source = Origin[key]
        elif STEP == 6:
            Link = self.Link_Memo_fastq
            Link_Memo = self.Link_Memo_gsnap
            with open(Link, "rb") as f1:
                Link_source = pickle.load(f1)

        with open(Link_Memo, "rb") as f:
            DESTINATION = pickle.load(f)
            for origin in Link_source:
                for destination in DESTINATION:
                    if (origin.split("/")[-1] == destination.split("/")[-1]) == True:
                        try:
                            print(origin, "is now linking with", destination)
                            os.symlink(origin, destination)
                        except FileExistsError as er:
                            print(er)
                            continue
                    else:
                        continue

    def Rename(self, A, B):
        os.system('mv %s %s' % (A, B))
        print("%s >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> \n %s " % (A, B))



if __name__ == "__main__":
    #exec_RNA("20170601", "RSq")

    trial = Rename_Center("20170919")
    query = "1330"
    replaceFROM = "_RSq"
    replaceTO = "_T01_RSq"

    ####### Step 1 ################isoform dir 3개만 먼저 파일명, 폴더명을 바꾼다 isoform dir의 gsnap.gz 링크를 해제한다.링크 파일 이름만 기억시켜둔다##############
    #trial.Rename_Subdirectories_EQL8(query, replaceFROM, replaceTO, "Files", 1)
    #trial.Rename_Subdirectories_EQL8(query, replaceFROM, replaceTO, "Dir", 1)
    #trial.making_Unlink(query, 1)

    ########Step 2########## mut dir의 파일명, 폴더명을 교체한다. fq.gz 링크를 해제한다.####링크 파일이름만 기억시켜둔다######################

    #trial.Rename_Subdirectories_EQL8(query, replaceFROM, replaceTO, "Files", 2)
    #trial.Rename_Subdirectories_EQL8(query, replaceFROM, replaceTO, "Dir", 2)
    #trial.making_Unlink(query, 2)

    ##########Step 3####### expr dir의 파일, 폴더명을 교체한다. fq.gz를 해제한다. 링크 파일 이름만 기억시켜둔다.
    #trial.Rename_Subdirectories_EQL8(query, replaceFROM, replaceTO, "Files", 3)
    #trial.Rename_Subdirectories_EQL8(query, replaceFROM, replaceTO, "Dir", 3)
    #trial.making_Unlink(query, 3)

    ###################### isoform dir의 gsnap.gz링크를 다시 걸어준다 ##########################
    #trial.construct_Link(query, 4)

    ################## EQL2의 링크 이름을 바꿔준다. expr dir, mut dir 순대로  fq.gz링크를 걸어준다
    #trial.Rename_Rootdirectories_EQL2(query, replaceFROM, replaceTO)
    #trial.construct_Link(query, 5)
    #trial.construct_Link(query, 6)


    ############################위의 과정을 축약한 것이 아래의 코드이다##################################################################

    ######## ALL STEPS AT ONCE ######
    for i in range(1,4,1):
        print("Step", i)
        trial.Rename_Subdirectories_EQL8(query, replaceFROM, replaceTO, "Files", i)
        trial.Rename_Subdirectories_EQL8(query, replaceFROM, replaceTO, "Dir", i)
        trial.making_Unlink(query, i)

    for j in range(4,7,1):
        print("Step", j)
        trial.construct_Link(query, j)
        if j == 4:
            trial.Rename_Rootdirectories_EQL2(query, replaceFROM, replaceTO)
        else:
            continue
