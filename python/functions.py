import sys; sys.path.append('/opt/root/lib')
import ROOT
from ROOT import vector
from pyspark import SparkFiles

class RootLoad :
    class __Libraries :
        def __init__(self) :
            print "Loading Root Classes..."
            ROOT.gSystem.Load("EnergyLoss_C.so")
            ROOT.gSystem.Load("Calibrations_C.so")
            ROOT.Calibrations.InitParameters()
            ROOT.gSystem.Load("EnergyAngle_C.so")
            ROOT.EnergyAngle.ReadLookupTable()
            ROOT.gSystem.Load("Spectra_C.so")
            from ROOT import GoodEvent
    instance_Libraries = None
    def __init__(self) :
        if not RootLoad.instance_Libraries :
            RootLoad.instance_Libraries = RootLoad.__Libraries()
    def __getattr__(self,name) :
        return getattr(self.instance_Libraries,name)

def process_energy_angle(index,iterat) :
    RootLoad()
    #Extract run number
    filename = iterat.next()
    run = int(filename[1][-13:-10])
    print "Opening run {0}".format(run)
    #open root file, make new instance of EnergyAngle class
    root_file = ROOT.THDFSFile(filename[1],"hdfs://cycdhcp22.tamu.edu:54310")
    root_tree = root_file.Get("rawData")
    reader = ROOT.EnergyAngle(root_tree)
    #loop tree, filling values in EnergyAngle
    reader.Loop(index)
    return [index]

def get_scattering_events(index) :
    RootLoad()
    file = ROOT.TFile("energy_angle_{0}.root".format(index))
    tree = file.Get("energyAngle")
    spectra = ROOT.Spectra(tree)
    good_events = spectra.Loop()
    good_event_list = list()
    for e in good_events :
        this_event = dict()
        this_event['pc_wire'] = e.wire
        this_event['si_label']="{0}-{1}".format(e.detector,e.quadrant)
        this_event['cm_energy'] = e.cm_energy
        this_event['pc_position'] = e.position
        good_event_list.append(this_event)
    return good_event_list
 
