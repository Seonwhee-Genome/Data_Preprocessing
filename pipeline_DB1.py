import json, pickle
import sqlite3
import pandas as pd
import pymysql
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
        jsons = glob("/home/jsgene/JSONS/CHA*")
        NGS = self.connect_mongoDB()

        #for ajson in jsons:
        #    A_json = json.load(open(ajson))
        #    BTS_ID = NGS.insert_one(A_json).inserted_id
        #    print(NGS.find_one({"_id": BTS_ID}))
        BTS_ID = NGS.insert_many(jsons)

    def slicing_DB(self):
        NGS = self.connect_mongoDB()
        for post in NGS.find({"ircr_name": "BT17_1353_B"}):
            print(post)







if __name__=="__main__":
    #MySQL_connection('ircr1', 'rpkm_subtype_ssgsea')
    #MySQL_connection('ircr1', 'mutation')
    #MySQL_connection('ircr1', 'sample_info')
    #MySQL_connection('ircr1', 'rpkm_gene_expr')
    nosql = NO_SQL()
    #nosql.connect_mongoDB()
    nosql.slicing_DB()
