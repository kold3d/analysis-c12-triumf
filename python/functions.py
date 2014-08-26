import sys; sys.path.append('/opt/root/lib')
import ROOT
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
    instance_Libraries = None
    def __init__(self) :
        if not RootLoad.instance_Libraries :
            RootLoad.instance_Libraries = RootLoad.__Libraries()
    def __getattr__(self,name) :
        return getattr(self.instance_Libraries,name)

def process_file(filename) :
    RootLoad()
    #Extract run number
    run = int(filename[1][-13:-10])
    print "Opening run {0}".format(run)
    #open root file, make new instance of EnergyAngle class
    root_file = ROOT.THDFSFile(filename[1],"hdfs://cycdhcp22.tamu.edu:54310")
    root_tree = root_file.Get("rawData")
    reader = ROOT.EnergyAngle(root_tree)
    #loop tree, filling values in EnergyAngle
    nentries = reader.fChain.GetEntriesFast()
    event_list = list()
    for jentry in xrange(0,nentries) :
        ientry = reader.LoadTree(jentry)
        if ientry < 0 : break
        reader.fChain.GetEntry(jentry)
        #declare event dictionary
        if jentry % 100000 == 0 : 
            print "Processed {0} events of {1}".format(jentry,nentries)
        event_dictionary = dict();
        if ord(reader.si_mul) > 0 :
            si_dict = dict()
            for i in xrange(0,ord(reader.si_mul)) :
                det = ord(reader.si_det[i])
                quad = ord(reader.si_quad[i])
                cal_value = ROOT.Calibrations.CalibrateSi(reader.si_ch_e[i],det,quad)
                if cal_value > 350 : 
                    si_dict["{0}-{1}".format(det,quad)] = cal_value            
            event_dictionary['si-e']=si_dict
        if ord(reader.pc_mul) > 0 :
            pc_dict = dict()
            for i in xrange(0,ord(reader.pc_mul)) :
                wire = ord(reader.pc_wire[i])
                pc_dict["{0}-l".format(wire)]= \
                  ROOT.Calibrations.MatchPCLeft(reader.pc_ch_left_e[i],wire,run)
                pc_dict["{0}-r".format(wire)]= \
                  ROOT.Calibrations.MatchPCRight(reader.pc_ch_right_e[i],wire,run)
            event_dictionary['pc-e']=pc_dict
        event_dictionary['rf-t'] = reader.rf_t
        event_dictionary['ic-e'] = reader.ic_ch_e
        if len(event_dictionary) > 0 : event_list.append(event_dictionary)
    root_file.Close()
    print "There are {0} entries in the event list.".format(len(event_list))
    return event_list

def process_raw(raw) :
    event = dict()
    si_data = raw.get("si-e")
    pc_data = raw.get("pc-e")
    if si_data == None or len(si_data) == 0 or \
       pc_data == None or len(pc_data) == 0 : return event
    high_si_energy = 0.
    high_si_label = ''
    for key,value in si_data.iteritems() :
        if value > high_si_energy :
            high_si_energy = value
            high_si_label = key
    if high_si_label == '' : return event
    high_pc_energy = 0.
    high_pc_position = 0.
    high_pc_wire = 0
    for wire in xrange(1,9) :
        left_label = "{0}-l".format(wire)
        right_label = "{0}-r".format(wire)
        if pc_data.get(left_label) != None and pc_data.get(right_label) != None :
            pc_sum = pc_data[left_label]+pc_data[right_label]
            if pc_sum > high_pc_energy :
                high_pc_energy = pc_sum
                high_pc_wire = wire
                high_pc_position = ROOT.Calibrations.CalcPosition(wire-1,pc_data[left_label],pc_data[right_label])
    if high_pc_wire == 0 : return event
    event['si_energy']=high_si_energy
    event['si_label']=high_si_label
    event['pc_energy'] = high_pc_energy
    event['pc_wire'] = high_pc_wire
    event['pc_position'] = high_pc_position
    event['rf-t'] = raw['rf-t']
    event['ic-e'] = raw['ic-e']
    return event

def is_proton(event,cuts) :
    pc_energy = event.get("pc_energy")
    si_energy = event.get("si_energy")
    if pc_energy == None or si_energy == None :
        return False
    else :
        pc_wire = event.get("pc_wire")
        return cuts[pc_wire-1].IsInside(si_energy,pc_energy) and \
            cuts[8].IsInside(si_energy,event['rf-t']) and \
            event['ic-e'] < 230 and event['ic-e'] > 120

def lookup_cm_energy(event)  :
    lookup_result = ROOT.EnergyAngle.LookupCMEnergyAngle(event['pc_wire']-1,event['pc_position'],event['si_energy']/1000.)
    event['cm_energy'] = lookup_result.first
    event['cm_angle'] = lookup_result.second
    return event
