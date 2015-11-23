
from mut_framework.LJMET_converter.samples_25ns import data_samples,mc_samples
from mut_framework.BTagSF_FtCM.JetCounting import JetCountingManager
from ROOT import TFile
import  os
import argparse

taggers = {"pfCombinedInclusiveSecondaryVertexV2BJetTags": [0.604, 0.890, 0.970]}

default_input_dir = "/user_data/pdecastr/BTagSF/mut_ttrees/13TeV_25ns/selected/" 
default_output_dir = "../output/13TeV_25ns/" 
default_pt_bin_edges = [30, 50, 70, 100, 140]
default_eta_bin_edges = [-2.4, 2.4]

parser = argparse.ArgumentParser(description='integers.')
parser.add_argument('--input_folder', default=default_input_dir, type=str)
parser.add_argument('--output_folder', default=default_output_dir, type=str)
parser.add_argument('--mc_samples', default=[], nargs='*')
parser.add_argument('--data_samples', default=[], nargs='*')
parser.add_argument('--all_samples', action="store_true")
parser.add_argument('--pt_dependent', action="store_true")
parser.add_argument('--pt_bin_edges', default=default_pt_bin_edges, nargs='*')
parser.add_argument('--eta_bin_edges', default=default_eta_bin_edges, nargs='*')
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

if args.pt_dependent:    
    ptBins = args.pt_bin_edges
    etaBins = args.eta_bin_edges
else:
    ptBins = [0.,100000.]
    etaBins = [-2.4, 2.4]


mc_files = [ (i_folder+mc_samples[s]["full_name"]+".root") for s in mc_names]
data_files = [ (i_folder+data_samples[s]["full_name"]+".root") for s in data_names]

jetCM = JetCountingManager(ptBins, etaBins, taggers, oldMuonSF = False, out_dir=o_folder)  

jetCM.process(mc_files, isData = False )
jetCM.process(data_files, isData = True)
