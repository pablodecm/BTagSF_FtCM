
from mut_framework.BTagSF_FtCM.FitAssistant import FitAssistant
from ROOT import RooFit

fa = FitAssistant("../output/reduced_results.json")
fa.arg_list.Print()
fa.fix_parameters("combinedMVABJetTags", 0.4)
data = fa.data_hist("combinedMVABJetTags", 0.4)
fa.simulpdf.fitTo(data, RooFit.Extended())

