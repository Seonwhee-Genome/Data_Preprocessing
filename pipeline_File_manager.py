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
        print("Copying %s >>>>> To %s" %(FROM, TO))
        #sp.check_output(["cp", "-r", FROM, TO], stderr=sp.STDOUT, shell=True)
        sp.call("cp -r %s %s"%(FROM, TO), shell=True)

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

    def File_transfer(self):
        Old_location = "/EQL7/pipeline/GS_test"
        New_location = "/EQL7/pipeline/BS_test"
        Keywords = [ "S1318157_T", "S1318157_T", "IRCR_BT17_1319_B", "IRCR_BT15_870_B", "IRCR_BT15_870_T", "IRCR_BT16_1023_B", "IRCR_BT16_1023_T",
     "IRCR_GBM13_300_B", "IRCR_GBM13_300_T", "IRCR_GBM14_509_B", "IRCR_GBM14_509_T", "IRCR_GBM14_510_B",
     "IRCR_GBM14_510_T", "IRCR_BT16_935_B", "IRCR_BT16_935_T", "IRCR_BT15_925_B", "IRCR_BT15_925_T", "IRCR_BT17_1323_B",
     "IRCR_BT17_1323_T01", "IRCR_BT17_1323_T02", "IRCR_BT17_1323_T04", "IRCR_BT17_1325_B", "IRCR_BT17_1325_T01",
     "IRCR_BT17_1325_T03", "IRCR_GBM10_038_B", "IRCR_GBM10_038_T", "IRCR_GBM13_231_B", "IRCR_GBM13_231_T",
     "IRCR_GBM13_292_B", "IRCR_GBM13_292_T", "IRCR_GBM13_319_B", "IRCR_GBM13_319_T", "IRCR_GBM15_711_B",
     "IRCR_GBM15_711_T", "IRCR_GBM14_458_B", "IRCR_GBM14_458_T", "IRCR_GBM14_487_B", "IRCR_GBM14_487_T",
     "IRCR_GBM14_494_B", "IRCR_GBM14_494_T", "IRCR_GBM14_497_B", "IRCR_GBM14_497_T", "IRCR_GBM14_508_B",
     "IRCR_GBM14_508_T", "IRCR_GBM14_511_B", "IRCR_GBM14_511_T", "IRCR_GBM14_517_B", "IRCR_GBM14_517_T02",
     "IRCR_GBM14_571_B", "IRCR_GBM14_571_T", "IRCR_GBM14_614_B", "IRCR_GBM14_614_T", "IRCR_GBM14_616_B",
     "IRCR_GBM14_616_T", "IRCR_GBM14_618_B", "IRCR_GBM14_618_T", "IRCR_GBM15_686_B", "IRCR_GBM15_686_T",
     "IRCR_GBM14_534_B", "IRCR_GBM14_534_T", "IRCR_GBM14_589_B", "IRCR_GBM14_589_T", "IRCR_GBM14_567_B",
     "IRCR_GBM14_567_T", "IRCR_GBM14_586_B", "IRCR_GBM14_586_T", "IRCR_GBM14_626_B", "IRCR_GBM14_626_T",
     "CHA_BT17_023_B", "CHA_BT17_023_T03", "CHA_BT17_023_T04", "SNU_BT17_026_B", "SNU_BT17_026_T01", "SNU_BT17_027_B",
     "SNU_BT17_027_T", "IRCR_BT17_1379_B", "IRCR_BT17_1379_T", "IRCR_BT17_1384_B", "IRCR_BT17_1384_T",
     "IRCR_BT17_1351_B", "S16124067_T", "SNU_BT17_024_B", "SNU_BT17_024_T01", "SNU_BT17_025_B", "SNU_BT17_025_T01",
     "IRCR_BT17_1366_B", "IRCR_BT17_1366_T", "IRCR_BT17_1367_B", "IRCR_BT17_1367_T", "IRCR_BT17_1373_B",
     "IRCR_BT17_1373_T01", "IRCR_BT17_1373_T02", "IRCR_BT17_1373_T04", "IRCR_BT17_1375_B", "IRCR_BT17_1375_T",
     "IRCR_BT17_1340_B", "IRCR_BT17_1340_T", "SNU_BT17_028_B", "SNU_BT17_028_T", "CHA_BT17_029_B", "CHA_BT17_029_T",
     "IRCR_BT17_1389_B", "IRCR_BT17_1389_T01", "IRCR_BT17_1389_T02", "IRCR_BT17_1394_B", "IRCR_BT17_1394_T",
     "IRCR_BT17_1395_B", "IRCR_BT17_1395_T", "IRCR_BT17_1397_B", "IRCR_BT17_1397_T"]
        for word in Keywords:
            old_dir = "%s/%s_GS" %(Old_location, word)
            old_html = "%s.html" % old_dir
            new_dir = "%s/%s_BS" %(New_location, word)
            new_html = "%s.html" % new_dir
            self.move_whole_dir(old_dir, new_dir)
            self.move_whole_dir(old_html, new_html)

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

    ###########################
    ##mgr = File_manager(glob("/EQL11/pipeline/*"))
    ########################## Move Files to IRCR storage server #######################################################
    # mgr.from_EQL_to_Storage("RSq")

    ############################ Create Multiple Symbolic Links ########################################################
    ##mgr.Linking_multiple_dirs()

    ############################# Reset Symbolic Links and Create again ################################################

    #mgr = File_manager()
    #mgr.File_transfer()
    #mgr.Link_modifier("TERAGEN", "C00", "20160826", "RSq)
    #mgr.Link_modifier("GS", "", "20170607", "GS")
