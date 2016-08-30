#python scripts/makejobs.py --outputFolder jobs --fileList filelists/ZeroBias_RunBAndE.dat --splitting 30 --outputFilename output/RateBinnedByLumi --options "--runs 273158,273302,273402,273403,273404,273405,273406,273408,273409,273410,273411,273425,277069,277070,277071,277072,277073,277076,277087,277094,277096,277112,277126,277127,277148 --isData"

#! /usr/bin/env python

import os
import stat
import sys

from optparse import OptionParser

parser = OptionParser()
parser.add_option('--outputFolder',dest="outputFolder",help='')
parser.add_option('--fileList',    dest="fileList",help='')
parser.add_option('--options',     dest="options", help='')
parser.add_option('--splitting',   dest="splitting", help='')
parser.add_option('--outputFilename',   dest="outputFilename", help='')
(options, args) = parser.parse_args()

print "outputFolder: ",options.outputFolder
print "fileList    : ",options.fileList
print "options     : ",options.options
print "splitting   : ",options.splitting
print "outputFilename: ",options.outputFilename

# If the output folder does not exist create it
if options.outputFolder is not None:
  outputFolder=options.outputFolder+"/"
  if not os.path.exists(outputFolder):
    os.makedirs(outputFolder)

import ROOT
chain = ROOT.TChain("l1EventTree/L1EventTree")

num_lines = sum(1 for line in open(options.fileList))


fSubmit = open(options.outputFolder+'/submitJobs.sh','w')
fSubmit.write("#!/bin/sh\n")
fSubmit.write("\n")

fHadd = open(options.outputFolder+'/haddJobs.sh','w')
fHadd.write("#!/bin/sh\n")
fHadd.write("hadd -f " + options.outputFilename + ".root")

import math

for x in range(0,int(math.ceil(num_lines/float(options.splitting)))):
  
  minFile = str(x * int(options.splitting))
  maxFile = str((x+1) * int(options.splitting))
  
  fSubmit.write("qsub -q hep.q -l h_rt=0:30:0 " + options.outputFolder + "/runJob_"+str(x)+".sh\n")
  
  fHadd.write(" " + options.outputFilename + str(x)+".root")
  
  with open(outputFolder+'runJob_'+str(x)+'.sh','w') as fOut:
    fOut.write("#!/bin/bash\n")
    fOut.write("\n")
    fOut.write("cd /home/hep/dw515/NTupleProduction/l1t-integration-v62.0/CMSSW_8_0_9/src/TriggerStudies/L1VBFTrigger/\n");
    fOut.write("source /vols/grid/cms/setup.sh\n");
    fOut.write("eval `scramv1 runtime -sh`\n");
    fOut.write("EventsPassedL1 --input "+options.fileList+" --outputFilename " + options.outputFilename + str(x)+".root "+" --minFile " + minFile + " --maxFile " + maxFile +" "+ options.options+"\n")

  os.chmod(outputFolder+'/runJob_'+str(x)+'.sh', stat.S_IRUSR | stat.S_IWUSR | stat.S_IRGRP | stat.S_IROTH | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)

fHadd.write("\n")
for x in range(0,int(math.ceil(num_lines/float(options.splitting)))):
  fHadd.write("rm " + options.outputFilename + str(x)+".root\n")
  
os.chmod(outputFolder+'/submitJobs.sh', stat.S_IRUSR | stat.S_IWUSR | stat.S_IRGRP | stat.S_IROTH | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)
os.chmod(outputFolder+'/haddJobs.sh', stat.S_IRUSR | stat.S_IWUSR | stat.S_IRGRP | stat.S_IROTH | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)
