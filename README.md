# Data_Preprocessing
Data calling, saving, filtering from either DB or local files

### Installation of mongoDB
[root@smc-centos jsgene]# mv mongodb-linux-x86_64-2.6.6.tgz /usr/local
[root@smc-centos jsgene]# cd /usr/local
[root@smc-centos local]# tar -zxvf mongodb-linux-x86_64-2.6.6.tgz  
mongodb-linux-x86_64-2.6.6/README
mongodb-linux-x86_64-2.6.6/THIRD-PARTY-NOTICES
mongodb-linux-x86_64-2.6.6/GNU-AGPL-3.0
mongodb-linux-x86_64-2.6.6/bin/mongodump
mongodb-linux-x86_64-2.6.6/bin/mongorestore
mongodb-linux-x86_64-2.6.6/bin/mongoexport
mongodb-linux-x86_64-2.6.6/bin/mongoimport
mongodb-linux-x86_64-2.6.6/bin/mongostat
mongodb-linux-x86_64-2.6.6/bin/mongotop
mongodb-linux-x86_64-2.6.6/bin/mongooplog
mongodb-linux-x86_64-2.6.6/bin/mongofiles
mongodb-linux-x86_64-2.6.6/bin/bsondump
mongodb-linux-x86_64-2.6.6/bin/mongoperf
mongodb-linux-x86_64-2.6.6/bin/mongod
mongodb-linux-x86_64-2.6.6/bin/mongos
mongodb-linux-x86_64-2.6.6/bin/mong

[root@smc-centos local]# mv mongodb-linux-x86_64-2.6.6 mongodb
[root@smc-centos local]# cd mongodb

[root@smc-centos mongodb]# mkdir data
[root@smc-centos mongodb]# mkdir config
[root@smc-centos mongodb]# mkdir log
[root@smc-centos mongodb]# cd config/

[root@smc-centos config]# vi mongodb.conf
dbpath=/usr/local/mongodb
logpath=/usr/local/mongodb/log/mongodb.log
logappend=true
port=5050
verbose=true
fork=true
rest=true

[root@smc-centos config]# cd ..
[root@smc-centos mongodb]# cd bin
[root@smc-centos bin]# ./mongod --config /usr/local/mongodb/config/mongodb.conf 

2017-08-30T12:51:05.817+0900 ** WARNING: --rest is specified without --httpinterface,
2017-08-30T12:51:05.817+0900 **          enabling http interface
about to fork child process, waiting until server is ready for connections.
forked process: 20010
child process started successfully, parent exiting
 
[root@smc-centos bin]# ./mongo localhost:5050
MongoDB shell version: 2.6.6
connecting to: localhost:5050/test
Welcome to the MongoDB shell.
For interactive help, type "help".
For more comprehensive documentation, see
        http://docs.mongodb.org/
Questions? Try the support group
        http://groups.google.com/group/mongodb-user
Server has startup warnings: 
2017-08-30T12:51:05.817+0900 ** WARNING: --rest is specified without --httpinterface,
2017-08-30T12:51:05.817+0900 **          enabling http interface
2017-08-30T12:51:05.834+0900 [initandlisten] 
2017-08-30T12:51:05.835+0900 [initandlisten] ** WARNING: You are running on a NUMA machine.
2017-08-30T12:51:05.835+0900 [initandlisten] **          We suggest launching mongod like this to avoid performance problems:
2017-08-30T12:51:05.835+0900 [initandlisten] **              numactl --interleave=all mongod [other options]
2017-08-30T12:51:05.835+0900 [initandlisten] 
> 
