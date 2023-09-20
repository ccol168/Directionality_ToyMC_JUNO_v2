import os
import random
import sys
import argparse

def ManageWarning () :
    while True :
        data = input()
        if data == "yes" or data == "y" :
            return
        elif data == "no" or data == "n" :
            exit(1)
        else : 
            print("Invalid response, type yes or no")
            continue

#PARSING
prs = argparse.ArgumentParser()           # parser name

prs.add_argument("-name","--NameRun",help="Name of the run",required=True)
prs.add_argument("-runs", "--Runs", help="How many processes to launch",type=int,required=True)
prs.add_argument("-events", "--EventsPerRun", help="How many events per run",type=int,required=True)
prs.add_argument("-cfg", "--CfgFile", help="Path to the configuration file template",required=True)


args = prs.parse_args()

NameRun = str(args.NameRun)
Runs = args.Runs
EventsPerRun = args.EventsPerRun
CfgFile=str(args.CfgFile)

if os.path.exists("Events/" + NameRun) :
    print("WARNING : A run named " + NameRun + " already exists, want to continue anyway?")
    ManageWarning()

if not os.path.exists("Events") :
    print("ERROR: the program needs an Event directory to store results")
    exit(1)

if not os.path.exists(CfgFile) :
    print("ERROR: Configuration file not found")
    exit(1)



#Creating Directories
if not os.path.exists("Events/" + NameRun):
    os.makedirs("Events/" + NameRun)

if not os.path.exists("Events/" + NameRun + '/log'):
    os.makedirs("Events/" + NameRun + '/log')

if not os.path.exists("Events/" + NameRun + '/err'):
    os.makedirs("Events/" + NameRun + '/err')

if not os.path.exists("Events/" + NameRun + '/out'):
    os.makedirs("Events/" + NameRun + '/out')

if not os.path.exists("Events/" + NameRun + '/root'):
    os.makedirs("Events/" + NameRun + '/root')

if not os.path.exists("Events/" + NameRun + '/sh'):
    os.makedirs("Events/" + NameRun + '/sh')

if not os.path.exists("Events/" + NameRun + '/sub'):
    os.makedirs("Events/" + NameRun + '/sub')

if not os.path.exists("Events/" + NameRun + '/configs'):
    os.makedirs("Events/" + NameRun + '/configs')

if not os.path.exists("Events/" + NameRun + '/analyzed'):
    os.makedirs("Events/" + NameRun + '/analyzed')

#Reading the configuration file and copying it for each run
for i in range(0,Runs) :
    ReadCfgFile = open(CfgFile,'r')
    WriteFile = open("Events/" + NameRun + "/configs/CfgFile_" +str(i).zfill(4)+".txt","w")

    for j in range (0,2) :
        line = ReadCfgFile.readline()
        WriteFile.write(line)

    ReadCfgFile.readline()
    WriteFile.write("NEvents " + str(EventsPerRun) + "\n")

    for j in range(0,13) :
        line = ReadCfgFile.readline()
        WriteFile.write(line)

    #reading if it is backgrounds or not    
    line = ReadCfgFile.readline()
    WriteFile.write(line)
    content = line.split()
    IsBackgrounds = content[1]

    #reading the last cut
    line = ReadCfgFile.readline()
    WriteFile.write(line)

    ReadCfgFile.close()
    WriteFile.close()

DirPath = os.getcwd() + "/Events/" + NameRun + "/"

#writing the .sh files
for i in range(0,Runs) :
    WriteFile = open("Events/" + NameRun + "/sh/Script_" +str(i).zfill(4)+".sh","w")
    WriteFile.write('#!/bin/bash \n')
    WriteFile.write('source /cvmfs/juno.ihep.ac.cn/centos7_amd64_gcc1120/Pre-Release/J23.1.0-rc0/setup.sh \n \n')

    WriteFile.write(os.getcwd() + '/Directionality_ToyMC  ')
    WriteFile.write(DirPath + 'configs/CfgFile_' + str(i).zfill(4) + '.txt  ' )
    WriteFile.write(DirPath + 'root/' + NameRun + '_' + str(i).zfill(4) + '.root \n \n')

    WriteFile.write('sleep 10 \n \n')

    WriteFile.write(os.getcwd() + '/Ordered_hits  ')
    WriteFile.write(DirPath + 'root/' + NameRun + '_' + str(i).zfill(4) + '.root  ')
    WriteFile.write(DirPath + 'analyzed/Analyzed_' + NameRun + '_' + str(i).zfill(4) + '.root  60 ' + IsBackgrounds)

    WriteFile.close()

#writing the .sub files 

for i in range(0,Runs) :
    WriteFile = open("Events/" + NameRun + "/sub/Script_" +str(i).zfill(4)+".sub","w")
    WriteFile.write('universe = vanilla\n')
    WriteFile.write('getenv = true\n')
    WriteFile.write('executable = ' + os.getcwd() + "/Events/" + NameRun + "/sh/Script_" + str(i).zfill(4) + '.sh \n')
    WriteFile.write('log = ' + os.getcwd() + "/Events/" + NameRun + '/log/' + NameRun + '_' + str(i).zfill(4) + '.log\n')
    WriteFile.write('output = ' + os.getcwd() + "/Events/" + NameRun + '/out/' + NameRun + '_' + str(i).zfill(4)  + '.out\n')
    WriteFile.write('error = ' + os.getcwd() + "/Events/" + NameRun + '/err/' + NameRun + '_' + str(i).zfill(4) + '.err\n')
    WriteFile.write('+MaxRuntime = 86400\n')
    WriteFile.write('ShouldTransferFiles = YES\n')
    WriteFile.write('WhenToTransferOutput = ON_EXIT\n')
    WriteFile.write('queue 1\n')
    WriteFile.close()

#writing the CondorScriptToLaunch.sh 


WriteFile = open("Events/" + NameRun + '/CondorScriptToLaunch_' + NameRun + '.sh','w') 

for i in range(0,Runs) :
    WriteFile.write('condor_submit -name sn-02.cr.cnaf.infn.it ' + os.getcwd() + '/Events/' + NameRun + "/sub/Script_" +str(i).zfill(4)+".sub \n")
WriteFile.close()