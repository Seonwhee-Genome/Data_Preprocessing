import json, pickle
import sqlite3
import pandas as pd
import pymysql
import subprocess as sp
import os
from pymongo import MongoClient
from glob import glob

def MySQL_connection(DB_name, Table_name):

    conn = pymysql.connect(host='localhost', user='root', password='123456', db='%s'%DB_name, charset='utf8')
    curs = conn.cursor()
    #curs = conn.cursor(pymysql.cursors.DictCursor)        # return Rows as dict, not as tuple
    sql_column_name = "SHOW COLUMNS FROM %s" %Table_name
    sql = "SELECT * FROM %s" %Table_name
    curs.execute(sql_column_name)
    headers = curs.fetchall()
    columnList = []
    rowList=[]

    curs.execute(sql)
    rows = curs.fetchall()
    for tuple in headers:
        columnList.append(tuple[0])

    for tuple in rows:
        rowList.append(tuple)

    df = pd.DataFrame(rowList, columns=columnList)
    print(df)
    with open("/data1/home/jsgene/PICKLES/%s_%s.pkl" %(DB_name, Table_name), 'w') as f:
        pickle.dump(df, f)

    conn.close()


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


    nosql = NO_SQL()
    #nosql.insert_data()
    nosql.Sample_lookup()
