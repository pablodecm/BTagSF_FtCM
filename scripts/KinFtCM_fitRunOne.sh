#!/bin/sh

source /cvmfs/cms.cern.ch/cmsset_default.sh
eval `scramv1 runtime -sh`
python KinFtCM_fitRunOne.py $1 $2
