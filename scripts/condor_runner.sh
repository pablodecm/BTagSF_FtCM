#!/bin/sh

source /cvmfs/cms.cern.ch/cmsset_default.sh
eval `scramv1 runtime -sh`
echo "Current dir:"$PWD
echo "Comand : "python $1 $2 $3 $4 $5 
python $1 $2 $3 $4 $5
