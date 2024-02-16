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

DirPath = os.getcwd() + "/Events/" + NameRun + "/"

#Reading the configuration file and copying it for each run
for i in range(0,Runs) :
    ReadCfgFile = open(CfgFile,'r')
    WriteFileCfg = open("Events/" + NameRun + "/configs/CfgFile_" +str(i).zfill(4)+".txt","w")

    for line in ReadCfgFile:
        line_content = line.split()

        if line_content[0] == "NEvents" :
            output_line = "NEvents " + str(EventsPerRun) + "\n"
        else :
            output_line = line
            if line_content[0] == "IsBackgrounds" :
                IsBackgrounds = line_content[1]
            if line_content[0] == "Fastmode" :
                Fastmode = line_content[1]
            if line_content[0] == "EventLimit" :
                MaxNEvents = line_content[1]
        WriteFileCfg.write(output_line)

    ReadCfgFile.close()
    WriteFileCfg.close()

    if Fastmode == '1' :
        EffectiveNEvents = MaxNEvents
    elif Fastmode == '0' :
        EffectiveNEvents = 60

#writing the .sh files
    WriteFileSh = open("Events/" + NameRun + "/sh/Script_" +str(i).zfill(4)+".sh","w")
    WriteFileSh.write('#!/bin/bash \n')
    WriteFileSh.write('source /cvmfs/juno.ihep.ac.cn/centos7_amd64_gcc1120/Pre-Release/J23.1.0-rc0/setup.sh \n \n')

    WriteFileSh.write(os.getcwd() + '/Directionality_ToyMC  ')
    WriteFileSh.write(DirPath + 'configs/CfgFile_' + str(i).zfill(4) + '.txt  ' )
    WriteFileSh.write(DirPath + 'root/' + NameRun + '_' + str(i).zfill(4) + '.root \n \n')

    WriteFileSh.write('sleep 10 \n \n')

    WriteFileSh.write(os.getcwd() + '/Ordered_hits  ')
    WriteFileSh.write(DirPath + 'root/' + NameRun + '_' + str(i).zfill(4) + '.root  ')
    WriteFileSh.write(DirPath + 'analyzed/Analyzed_' + NameRun + '_' + str(i).zfill(4) + '.root ' + EffectiveNEvents + ' ' + IsBackgrounds)

    WriteFileSh.close()

#writing the .sub files 

    WriteFileSub = open("Events/" + NameRun + "/sub/Script_" +str(i).zfill(4)+".sub","w")
    WriteFileSub.write('universe = vanilla\n')
    WriteFileSub.write('getenv = true\n')
    WriteFileSub.write('executable = ' + os.getcwd() + "/Events/" + NameRun + "/sh/Script_" + str(i).zfill(4) + '.sh \n')
    WriteFileSub.write('log = ' + os.getcwd() + "/Events/" + NameRun + '/log/' + NameRun + '_' + str(i).zfill(4) + '.log\n')
    WriteFileSub.write('output = ' + os.getcwd() + "/Events/" + NameRun + '/out/' + NameRun + '_' + str(i).zfill(4)  + '.out\n')
    WriteFileSub.write('error = ' + os.getcwd() + "/Events/" + NameRun + '/err/' + NameRun + '_' + str(i).zfill(4) + '.err\n')
    WriteFileSub.write('+MaxRuntime = 86400\n')
    WriteFileSub.write('ShouldTransferFiles = YES\n')
    WriteFileSub.write('WhenToTransferOutput = ON_EXIT\n')
    WriteFileSub.write('queue 1\n')
    WriteFileSub.close()

#writing the CondorScriptToLaunch.sh 
WriteFile = open("Events/" + NameRun + '/CondorScriptToLaunch_' + NameRun + '.sh','w') 

for i in range(0,Runs) :
    WriteFile.write('condor_submit -name sn-02.cr.cnaf.infn.it ' + os.getcwd() + '/Events/' + NameRun + "/sub/Script_" +str(i).zfill(4)+".sub \n")
WriteFile.close()