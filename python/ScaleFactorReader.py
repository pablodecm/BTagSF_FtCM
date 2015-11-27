from ROOT import BTagCalibration, BTagCalibrationReader

class ScaleFactorReader:
    def __init__(self, tagger, csv_file, op, meas = "comb"): 
        self.bc = BTagCalibration(tagger, csv_file)
        self.op = op
        self.meas = meas
    def get_sf( self, sys, flav, eta, pt):       
        br = BTagCalibrationReader(self.op, sys) 
        br.load(self.bc, flav, self.meas)
        return br.eval(flav, eta, pt)
