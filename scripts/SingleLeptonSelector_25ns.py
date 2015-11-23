#!/usr/bin/env python 

from mut_framework.LJMET_converter.samples_25ns import data_samples,mc_samples
from ROOT import TChain, SingleLeptonSelector 
import  os
import argparse


default_input_dir = "/user_data/pdecastr/BTagSF/mut_ttrees/13TeV_25ns/from_ljmet/" 
default_output_dir = "/user_data/pdecastr/BTagSF/mut_ttrees/13TeV_25ns/selected/" 

parser = argparse.ArgumentParser(description='integers.')
parser.add_argument('--input_folder', default=default_input_dir, type=str)
parser.add_argument('--output_folder', default=default_output_dir, type=str)
parser.add_argument('--mc_samples', default=[], nargs='*')
parser.add_argument('--data_samples', default=[], nargs='*')
parser.add_argument('--all_samples', action="store_true")
args = parser.parse_args()

i_folder = args.input_folder
o_folder = args.output_folder

if not os.path.exists(o_folder):
        os.makedirs(o_folder)

if args.all_samples:
    mc_names = mc_samples.keys()
    data_names = data_samples.keys()
else:
    mc_names = args.mc_samples
    data_names = args.data_samples

for name in mc_names:
    selector = SingleLeptonSelector()
    full_name = mc_samples[name]["full_name"]
    outfile = o_folder+full_name+".root" 
    tchain = TChain("tree")
    wildcard_name = i_folder + full_name + "*.root"
    print "Processing all files with wildcard:  {} ".format(wildcard_name)
    n_files = tchain.Add(wildcard_name)
    print "Added {} files to TChain".format(n_files)
    tchain.Process(selector,"ofile="+outfile)

for name in data_names:
    selector = SingleLeptonSelector()
    full_name = data_samples[name]["full_name"]
    outfile = o_folder+full_name+".root" 
    tchain = TChain("tree")
    wildcard_name = i_folder + full_name + "*.root"
    print "Processing all files with wildcard:  {} ".format(wildcard_name)
    n_files = tchain.Add(wildcard_name)
    print "Added {} files to TChain".format(n_files)
    tchain.Process(selector,"ofile="+outfile)
