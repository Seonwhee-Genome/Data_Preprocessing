import json, pickle
import sqlite3
import pandas as pd
import pymysql
import subprocess as sp
import os
from git import Repo, exc

from pymongo import MongoClient
from glob import glob

class SQL_DB(object):

    def __init__(self, DB_name, Table_name):
        self.connection = pymysql.connect(host='localhost', user='root', password='123456', db='%s' % DB_name, charset='utf8')
        self.cursor = self.connection.cursor()
        # self.connection.cursor(pymysql.cursors.DictCursor)        # return Rows as dict, not as tuple
        self.Table = Table_name
        self.DF_path = "/data1/home/jsgene/PICKLES/%s_%s" % (DB_name, Table_name)

    def Table_lookup(self, column_name, conditional=''):

        ################### Define Header ################################

        Header = []
        if column_name == "*":
            sql_column_name = "SHOW COLUMNS FROM %s" % self.Table
            self.cursor.execute(sql_column_name)
            headers = self.cursor.fetchall()
            for tuple in headers:
                Header.append(tuple[0])

        else:
            Headers = column_name.split(",")
            for head in Headers:
                Header.append("TCGA-%s" %head)
            print(Header)

        ##################################################################

        if conditional == '':
            query = "SELECT %s FROM %s" % (column_name, self.Table)
        else:
            query = "SELECT %s FROM %s WHERE %s" % (column_name, self.Table, conditional)
        self.cursor.execute(query)
        rows = self.cursor.fetchall()
        rowList = []

        for tuple in rows:
            rowList.append(tuple)

        ################## Output file ###################################
        df = pd.DataFrame(rowList, columns=Header)
        print(df)
        with open(self.DF_path+".pkl", 'w') as f:
            pickle.dump(df, f)

    def TCGA_rpkm(self, sample_names):
        ################### Define Header ################################
        Header = []
        Samples = []
        for aSample in sample_names.split(","):
            Samples.append("TCGA-%s" % aSample)
        print(Samples)
        sql_column_name = "SHOW COLUMNS FROM %s" % self.Table
        self.cursor.execute(sql_column_name)
        headers = self.cursor.fetchall()
        for tuple in headers:
            Header.append(tuple[0])
        List_Series = []
        ##################################################################
        for sample in Samples:
            query = "SELECT * FROM %s WHERE samp_id = '%s'" % (self.Table, sample)
            self.cursor.execute(query)
            rows = self.cursor.fetchall()
            Genes = []
            RPKMs = []
            for tuple in rows:
                Genes.append(tuple[1])
                RPKMs.append(tuple[-1])
            sample_rpkm = pd.Series(data=RPKMs, index=Genes, name=sample)
            List_Series.append(sample_rpkm)

        ################## Output file ###################################
        df = pd.concat(List_Series, axis=1)
        df.to_csv(self.DF_path+".csv")
        #with open(self.DF_path+".pkl", 'w') as f:
        #    pickle.dump(df, f)


    def TCGA_tables(self, sample_names):
        ################### Define Header ################################
        Header = []
        Samples = []
        for aSample in sample_names.split(","):
            Samples.append("TCGA-%s" % aSample)
        print(Samples)
        sql_column_name = "SHOW COLUMNS FROM %s" % self.Table
        self.cursor.execute(sql_column_name)
        headers = self.cursor.fetchall()
        for tuple in headers:
            Header.append(tuple[0])
        rowList = []
        ##################################################################
        for sample in Samples:
            query = "SELECT * FROM %s WHERE samp_id = '%s'" % (self.Table, sample)
            self.cursor.execute(query)
            rows = self.cursor.fetchall()

            for tuple in rows:
                rowList.append(tuple)


        ################## Output file ###################################
        df = pd.DataFrame(rowList, columns=Header)
        print(df)
        df.to_csv(self.DF_path+".csv")



    def __del__(self):
        self.connection.close()


class Meta_TABLE(object):
    def __init__(self):
        self.json_repository = "/data1/home/jsgene/JSONS"

    def json_loading(self, SampleName):
        A_json = json.load(open('%s/%s.json' %(self.json_repository, SampleName)))
        return A_json

    def build_SQL_TABLE(self):
        Sample_Dict = self.json_loading("IRCR_BT16_1236_T")
        someItem = Sample_Dict.itervalues().next()
        Columns = list(someItem.keys())
        Meta_DB = sqlite3.connect("%s/Metadata.sqlite")
        query = "INSERT INTO medicoes (TIMESTAMP,{0}) VALUES (?{1})"
        query = query.format(",".join(Columns), ",?" * len(Columns))
        print(query)

        for timestamp, data in Sample_Dict.iteritems():
            keys = (timestamp,) + tuple(data[c] for c in Columns)
            c = Meta_DB.cursor()
            c.execute(query)
            c.close()


class Git_Trace(object):

    def __init__(self, RepoPath):
        self.RepoPath = RepoPath
        new_repository = Repo.init(RepoPath)
        new_repository.close()
        self.repo = Repo(RepoPath)
        self.Index = self.repo.index
        self.Head = self.repo.head
        self.Git = self.repo.git
        self.Commit = self.repo.commit   # obtain commits at the specified revision

        with open('%s/.gitignore' %(self.RepoPath), 'w') as gi:
            gi.writelines("*/*bam\n*bam\n*/*bai\n*bai\n*/*gz\n*/*gsnap\n*/*zip")

    def __del__(self):
        print("Getting rid of Object %s from Memory!" %(self.repo))
        self.repo.close()

    def Creating_Branch(self):
        self.Git.checkout('HEAD', b="jsgene_new_branch")

    def Git_add(self, AddList):
        for group in AddList:
            add = self.Git.add(group)
            print(add)

    def Git_commit(self, Message):
        try:
            cmt = self.Git.commit(m=Message)
            print(cmt)
            print(self.Git.status())
        except exc.GitCommandError as sterr:
            print("nothing added to commit!!")


    def Get_commits(self):
        from time import strftime, gmtime
        fifty_commits = list(self.repo.iter_commits('master', max_count=50))
        for a_commit in fifty_commits:
            print(a_commit)
            print(a_commit.diff())
            print(strftime("%a, %d %b %Y %H:%M", gmtime(a_commit.committed_date)))
        headcommit = self.Head.commit
        print(headcommit.diff())

    def Git_delete(self):
        sp.call(["rm", "-rf", "%s/.git" %(self.RepoPath)])
        try:
            print(self.Git.status())
        except exc.GitCommandError as sterr:
            print("the git repository is successfully removed!!!")


class NO_SQL(object):
    def connect_mongoDB(self):
        client = MongoClient('localhost', 5050)
        DB = client.ircr1
        NGS = DB.ircr
        return NGS

    def insert_data(self):
        jsons = glob("/home/jsgene/JSONS/IRCR_BT17_1397*")
        NGS = self.connect_mongoDB()

        for ajson in jsons:
            A_json = json.load(open(ajson))
            BTS_ID = NGS.insert_one(A_json).inserted_id
            print(NGS.find_one({"_id": BTS_ID}))
        #BTS_ID = NGS.insert_many(jsons)

    def delete_data(self):
        name_to_delete = ["IRCR_BT15_1134_B", "IRCR_BT16_1154_B"]
        NGS = self.connect_mongoDB()
        for data in name_to_delete:
            print (data)
            target = {"ircr_name": {"$regex": data, "$options": "i"}}
            print("We got %s data" %NGS.count(target))
            result = NGS.delete_one(target)
            print("We remove %s data" %result.deleted_count)
            print("We got %s data" % NGS.count(target))



    def Sample_lookup(self):
        first_val = raw_input("ircr_name is ")
        NGS = self.connect_mongoDB()
        one_patient = {}
        paTH = {}
        initial_query = "%s" %first_val
        print ("query is %s" %initial_query)
        i = 1
        for post in NGS.find({"ircr_name": {"$regex" : initial_query, "$options": "i"}}):
            one_patient[i] = post["ircr_name"]
            i = i + 1

        try:
            print ("Please select number indicating one of your sample to lookup")
            print (one_patient)
            second_val = input("the key is ")
            print ("Here is the output >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>")
            print (one_patient[second_val])
            for sample in NGS.find({"ircr_name": one_patient[second_val]}):
                if sample["Scan_Type"] == "Glioma_Scan":
                    paTH[5] = "/EQL7/pipeline/GS_single/%s_GS/" % sample["new_ircr_name"]
                    paTH[4] = "/EQL7/pipeline/GS_test/%s_GS/" % sample["new_ircr_name"]
                    paTH[3] = "/EQL2/gbmscan_pool/%s_GS/" % sample["new_ircr_name"]
                    paTH[2] = "/EQL7/pipeline/%s/%s_GS/" % (sample["Group"], sample["new_ircr_name"])
                    paTH[1] = "/EQL2/%s/WXS/fastq/link/" % (sample["Group"])
                elif sample["Scan_Type"] == "BT_Scan":
                    paTH[5] = "/EQL7/pipeline/BS_single/%s_BS/" % sample["new_ircr_name"]
                    paTH[4] = "/EQL7/pipeline/BS_test/%s_BS/" % sample["new_ircr_name"]
                    paTH[3] = "/EQL2/btscan_pool/%s_BS/" % sample["new_ircr_name"]
                    paTH[2] = "/EQL7/pipeline/%s/%s_BS/" % (sample["Group"], sample["new_ircr_name"])
                    paTH[1] = "/EQL2/%s/WXS/fastq/link/" % (sample["Group"])

                elif sample["type"] == "RNA":
                    paTH[2] = "/EQL8/pipeline/%s/%s_RSq/" % (sample["Group"], sample["new_ircr_name"])
                    paTH[1] = "/EQL2/%s/RNASeq/fastq/link" % (sample["Group"])

            third_val = input("If you want to look up \n FASTQ files, please select 1. \n BAM files, please select 2 \n Copynumbers, please select 3 \n Variants, please select 4 \n")
            print ("Here is the output >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>")
            if third_val == 4:
                try:
                    sp.check_call(["ls %s" % paTH[third_val], "-al"], shell=True)
                except sp.CalledProcessError as err:
                    third_val = third_val+1  # =5
                    sp.check_call(["ls %s" % paTH[third_val], "-al"], shell=True)
            else:
                sp.check_call(["ls %s" % paTH[third_val], "-al"], shell=True)
            fileList = glob(paTH[third_val]+'*')
            j=1
            fileDict = {}
            for afile in fileList:
                fileDict[j] = afile
                j = j+1
            print(fileDict)
            forth_val = input("please select one file you want to see ")
            print ("Here is the output >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>")
            if fileDict[forth_val][-3:] == 'bam':
                (chrNum, chrStart, chrEnd) = tuple(int(x.strip()) for x in raw_input("You chose bam file, please type chromosome number and positions with comma seperation: (e.g : chr7:100677977-100678077 --> 7, 100677977, 100678077)").split(','))
                chrNum = "chr%s" %chrNum

                import pysam
                with pysam.Samfile(fileDict[forth_val], "rb") as bamfile:
                    for read in bamfile.fetch(chrNum, chrStart, chrEnd):
                        print (read)
            elif fileDict[forth_val][-6:] == '.fq.gz':
                import gzip
                from Bio import SeqIO
                Fastq = gzip.open(fileDict[forth_val], 'rb')
                print("We are now parsing a large Fastq file. Please wait")
                for record in SeqIO.parse(Fastq, "fastq"):
                    print(record.id)
                    print(record.seq)
                Fastq.close()

            else:
                sp.check_call(["cat %s" % fileDict[forth_val]], shell=True)

        except KeyError as ke:
            print (ke)
            print ("%s has not been registered yet!!" %val)



if __name__=="__main__":
    #MySQL_connection('ircr1', 'rpkm_subtype_ssgsea')
    #MySQL_connection('ircr1', 'mutation')
    #MySQL_connection('ircr1', 'sample_info')
    #MySQL_connection('ircr1', 'rpkm_gene_expr')
    GS_gene_set = ['SLC6A3', 'ERRFI1', 'SOX13', 'RB1', 'PLK1', 'NRAS', 'STK11', 'AKT1', 'AKT2', 'AKT3', 'ARAF', 'NOTCH4', 'NUP210L', 'MDM2', 'AXL', 'CREBZF', 'ZNF709', 'RBPJ', 'HGF', 'TRERF1', 'TTC30B', 'MDM4', 'PIK3CB', 'CACNA1C', 'PIK3CD', 'BCLAF1', 'FGR', 'PARK2', 'SEMG1', 'CDH8', 'SLC12A7', 'CBFB', 'MKLN1', 'CDH1', 'MYC', 'CACNA1S', 'ROS1', 'TOP2A', 'GLT8D2', 'BRCA1', 'CDKN2C', 'IGF1R', 'BRCA2', 'MED16', 'ERBB4', 'WRN', 'IDH2', 'IDH1', 'CXCR4', 'SAMD9L', 'ENTPD8', 'PIK3C2B', 'PIK3C2A', 'KEL', 'PIK3C2G', 'ARRDC1', 'ADAMTS20', 'FGFR4', 'SETD2', 'FGFR2', 'FGFR3','FGFR1', 'ZNF677', 'RSPO2', 'OPCML', 'ZC3H11A', 'SYK', 'ARID1B', 'CCND3', 'CCND2', 'FAM47C', 'ATF7IP2', 'GSK3B', 'ECHS1', 'GRB2', 'TLR6', 'MAX', 'APC', 'ZBTB20', 'NR1H3', 'LZTR1', 'ARID1A', 'PTPN11', 'G6PC', 'CD1D', 'PDGFA', 'CDKN1B', 'HRAS', 'FUBP1', 'MX2', 'PDGFB', 'CD209', 'PRG4', 'PFKFB1', 'SRPX', 'RERE', 'CD99L2', 'MET', 'GABRA6', 'EPHB4', 'DDX5', 'CPA2', 'THRA', 'TNRC18', 'CCNE1', 'STAT5A', 'SMO', 'ABHD13', 'PPAT', 'JAK2', 'H3F3B', 'NF1', 'ITK', 'BMX', 'ZNF766', 'CHEK2', 'JAK1', 'TGFBR2', 'CHEK1', 'IGF2', 'AFM', 'IL4R', 'CSK', 'ODF4', 'FAF1', 'ARID2', 'SENP5', 'TCF12', 'RNF43', 'PARP2', 'TNFRSF9', 'SLC26A3', 'VHL', 'MYO16', 'MUC17', 'HCK', 'FGF1', 'TYRP1', 'CCDC135', 'TP53', 'GNAQ', 'ESR1', 'AZGP1', 'MPL', 'DDR1', 'ATAD1', 'TOP1', 'KRT15', 'KRT13', 'TMEM216', 'PHF8', 'PIK3CA', 'EFCAB4B', 'PCDH19', 'FGF23', 'BTK', 'ZNF292', 'PDGFRB', 'PDGFRA', 'MAP3K1', 'GIGYF2', 'HDAC9', 'PDHA1', 'ATM', 'EPHA3', 'KRAS', 'RET', 'NTRK1', 'PTEN', 'CSMD3', 'HTRA2', 'ATR', 'MAP3K4', 'ABL1', 'FLT4', 'ABL2', 'RAF1', 'FLT3', 'BCR', 'FLT1', 'TEK', 'STAT3', 'STAT1', 'NPM1', 'MYH8', 'PTK6', 'CSF1R', 'TPTE2', 'PIK3R2', 'MAPK1', 'ZNF512B', 'MAPK8', 'DOT1L', 'DLC1', 'KDR', 'FBXW7', 'SRC', 'ETFB', 'MACC1', 'WWC2', 'HDLBP', 'CHRDL2', 'OR8B2', 'FOXL2', 'BCOR', 'MOB2', 'ATRX', 'JOSD2', 'ROCK2', 'ZDHHC4', 'PORCN', 'IFNA10', 'ALK', 'SMAD4', 'IRS4', 'ETV6', 'PLCG1', 'DDR2', 'SMC1A', 'EZH2', 'CDK4', 'CDK6', 'RAD50', 'HSP90AA1', 'ZMIZ1', 'PAN3', 'GNAS', 'CIC', 'AOX1', 'RNF168','RPL5', 'ERBB3', 'TRPV6', 'EEF1A1', 'CHGB', 'DIP2C', 'EGFR', 'HNF1A', 'BAX', 'CTNNB1', 'MAP2K2', 'MAP2K1', 'CDKN2A', 'FIZ1', 'MCL1', 'PIK3CG', 'MAP2K4', 'SOX2', 'NTRK2', 'SMARCA4', 'FAM83D', 'JAK3', 'QKI', 'RICTOR', 'PIK3R1', 'NAP1L2', 'RHBG', 'NPAS3', 'MLH1', 'FAM155A', 'KLF6', 'NIPBL', 'CTU1', 'AURKA', 'AURKB', 'AURKC', 'SIGLEC8', 'CNOT1', 'TERT', 'WNK3', 'DPP6', 'STK19', 'UTS2', 'NF2', 'KMT2A', 'DNMT3A', 'SMARCB1', 'PRDM2', 'MYCN', 'PRKCA', 'CD44', 'KANK1', 'PRKCD', 'PTCH2', 'AR', 'TSC2', 'LRRC55', 'CDKN2B', 'FAM53B', 'ST6GALNAC2', 'SPTA1', 'RNLS', 'THAP4', 'ROBO3', 'REN', 'ERBB2', 'CALCR', 'KIT', 'NOTCH1', 'EZR', 'NOTCH3', 'ZNF544', 'NOTCH2', 'BRAF', 'EMG1', 'GNA11', 'PASD1', 'MTOR', 'BAP1', 'PTCH1', 'ZNF548', 'BCL2']
    CS_gene_set = ['SLC6A3', 'ERRFI1', 'SOX13', 'RB1', 'PLK1', 'NRAS', 'STK11', 'AKT1', 'AKT2', 'AKT3', 'ARAF', 'NOTCH4', 'NUP210L', 'MDM2', 'AXL', 'CREBZF', 'ZNF709', 'RBPJ', 'HGF', 'TRERF1', 'TTC30B', 'MDM4', 'PIK3CB', 'CACNA1C', 'PIK3CD', 'BCLAF1', 'FGR', 'PARK2', 'SEMG1', 'CDH8', 'SLC12A7', 'CBFB', 'MKLN1', 'CDH1', 'MYC', 'CACNA1S', 'ROS1', 'TOP2A', 'GLT8D2', 'BRCA1', 'CDKN2C', 'IGF1R', 'BRCA2', 'MED16', 'ERBB4', 'WRN', 'IDH2', 'IDH1', 'CXCR4', 'SAMD9L', 'ENTPD8', 'PIK3C2B', 'PIK3C2A', 'KEL', 'PIK3C2G', 'ARRDC1', 'ADAMTS20', 'FGFR4', 'SETD2', 'FGFR2', 'FGFR3', 'FGFR1', 'ZNF677', 'RSPO2', 'OPCML', 'ZC3H11A', 'SYK', 'ARID1B', 'CCND3', 'CCND2', 'FAM47C', 'ATF7IP2', 'GSK3B', 'ECHS1', 'GRB2', 'TLR6', 'MAX', 'APC', 'ZBTB20', 'NR1H3', 'LZTR1', 'ARID1A', 'PTPN11', 'G6PC', 'CD1D', 'PDGFA', 'CDKN1B', 'HRAS', 'FUBP1', 'MX2', 'PDGFB', 'CD209', 'PRG4', 'PFKFB1', 'SRPX', 'RERE', 'CD99L2', 'MET', 'GABRA6', 'EPHB4', 'DDX5', 'CPA2', 'THRA', 'TNRC18', 'CCNE1', 'STAT5A', 'SMO', 'ABHD13', 'PPAT', 'JAK2', 'H3F3B', 'NF1', 'ITK', 'BMX', 'ZNF766', 'CHEK2', 'JAK1', 'TGFBR2', 'CHEK1', 'IGF2', 'AFM', 'IL4R', 'CSK', 'ODF4', 'FAF1', 'ARID2', 'SENP5', 'TCF12', 'RNF43', 'PARP2', 'TNFRSF9', 'SLC26A3', 'VHL', 'MYO16', 'MUC17', 'HCK', 'FGF1', 'TYRP1', 'CCDC135', 'TP53', 'GNAQ', 'ESR1', 'AZGP1', 'MPL', 'DDR1', 'ATAD1', 'TOP1', 'KRT15', 'KRT13', 'TMEM216', 'PHF8', 'PIK3CA', 'EFCAB4B', 'PCDH19', 'FGF23', 'BTK', 'ZNF292', 'PDGFRB', 'PDGFRA', 'MAP3K1', 'GIGYF2', 'HDAC9', 'PDHA1', 'ATM', 'EPHA3', 'KRAS', 'RET', 'NTRK1', 'PTEN', 'CSMD3', 'HTRA2', 'ATR', 'MAP3K4', 'ABL1', 'FLT4', 'ABL2', 'RAF1', 'FLT3', 'BCR', 'FLT1', 'TEK', 'STAT3', 'STAT1', 'NPM1', 'MYH8', 'PTK6', 'CSF1R', 'TPTE2', 'PIK3R2', 'MAPK1', 'ZNF512B', 'MAPK8', 'DOT1L', 'DLC1', 'KDR', 'FBXW7', 'SRC', 'ETFB', 'MACC1', 'WWC2', 'HDLBP', 'CHRDL2', 'OR8B2', 'FOXL2', 'BCOR', 'MOB2', 'ATRX', 'JOSD2', 'ROCK2', 'ZDHHC4', 'PORCN', 'IFNA10', 'ALK', 'SMAD4', 'IRS4', 'ETV6', 'PLCG1', 'DDR2', 'SMC1A', 'EZH2', 'CDK4', 'CDK6', 'RAD50', 'HSP90AA1', 'ZMIZ1', 'PAN3', 'GNAS', 'CIC', 'AOX1', 'RNF168', 'RPL5', 'ERBB3', 'TRPV6', 'EEF1A1', 'CHGB', 'DIP2C', 'EGFR', 'HNF1A', 'BAX', 'CTNNB1', 'MAP2K2', 'MAP2K1', 'CDKN2A', 'FIZ1', 'MCL1', 'PIK3CG', 'MAP2K4', 'SOX2', 'NTRK2', 'SMARCA4', 'FAM83D', 'JAK3', 'QKI', 'RICTOR', 'PIK3R1', 'NAP1L2', 'RHBG', 'NPAS3', 'MLH1', 'FAM155A', 'KLF6', 'NIPBL', 'CTU1', 'AURKA', 'AURKB', 'AURKC', 'SIGLEC8', 'CNOT1', 'TERT', 'WNK3', 'DPP6', 'STK19', 'UTS2', 'NF2', 'KMT2A', 'DNMT3A', 'SMARCB1', 'PRDM2', 'MYCN', 'PRKCA', 'CD44', 'KANK1', 'PRKCD', 'PTCH2', 'AR', 'TSC2', 'LRRC55', 'CDKN2B', 'FAM53B', 'ST6GALNAC2', 'SPTA1', 'RNLS', 'THAP4', 'ROBO3', 'REN', 'ERBB2', 'CALCR', 'KIT', 'NOTCH1', 'EZR', 'NOTCH3', 'ZNF544', 'NOTCH2', 'BRAF', 'EMG1', 'GNA11', 'PASD1', 'MTOR', 'BAP1', 'PTCH1', 'ZNF548', 'BCL2']


   
    #nosql = NO_SQL()
    #nosql.insert_data()
    #nosql.delete_data()
    #nosql.Sample_lookup()

    ###gits = Git_Trace('/data1/home/jsgene/TEST')
    ###gits.Git_add(["*txt"])
    ###gits.Git_commit("First Commit!!!")
    ###gits.Get_commits()
    #gits.Git_delete()

    sq = SQL_DB("tcga1", "splice_normal")
    #sq.Table_lookup("*", conditional="samp_id = 'TCGA-27-1830'")
    #sq.TCGA_rpkm("02-0047,06-0157,06-0158,06-0168,06-0171,06-0174,06-0178,06-0184,06-0187,06-0190,06-0210,06-0221,06-0238,06-0644,06-0645,06-0646,06-0649,06-0878,06-2570,06-5408,06-5410,06-5412,06-5413,06-5417,12-0616,12-3650,14-0736,14-0789,14-0790,14-0817,14-0871,14-1034,14-1402,14-1823,14-1825,14-1829,19-0957,19-1389,19-1390,19-1787,19-2619,19-2620,19-2624,19-2625,19-2629,19-4065,19-5960,27-1830,27-1831,27-1834,27-1835,27-1837,27-2519,27-2521,27-2523,27-2524,27-2526,27-2528,76-4925,76-4926,76-4927,76-4928,76-4929,76-4931,76-4932")
    sq.TCGA_tables("02-0047,06-0157,06-0158,06-0168,06-0171,06-0174,06-0178,06-0184,06-0187,06-0190,06-0210,06-0221,06-0238,06-0644,06-0645,06-0646,06-0649,06-0878,06-2570,06-5408,06-5410,06-5412,06-5413,06-5417,12-0616,12-3650,14-0736,14-0789,14-0790,14-0817,14-0871,14-1034,14-1402,14-1823,14-1825,14-1829,19-0957,19-1389,19-1390,19-1787,19-2619,19-2620,19-2624,19-2625,19-2629,19-4065,19-5960,27-1830,27-1831,27-1834,27-1835,27-1837,27-2519,27-2521,27-2523,27-2524,27-2526,27-2528,76-4925,76-4926,76-4927,76-4928,76-4929,76-4931,76-4932")




