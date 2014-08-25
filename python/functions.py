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

def process_si(data_string) :
    detector_dictionary = dict()
    events = data_string.split(',')
    for event in events :
        label_value_pair = event.split(':')
        ch_label = label_value_pair[0]
        ch_value = float(label_value_pair[1])
	cal_value = ROOT.Calibrations.CalibrateSi(ch_value,int(ch_label[0])-1,int(ch_label[2])-1)
        if cal_value > 350. : detector_dictionary[ch_label]=cal_value      
    return detector_dictionary
    
def process_pc(data_string,run) :
    detector_dictionary = dict()
    events = data_string.split(',')
    for event in events :
        label_value_pair = event.split(':')
        ch_label = label_value_pair[0]
        ch_value = float(label_value_pair[1])
        cal_value = 0.
        if ch_label[2] == 'l' : cal_value = ROOT.Calibrations.MatchPCLeft(ch_value,int(ch_label[0])-1,run)
        elif ch_label[2] == 'r' : cal_value = ROOT.Calibrations.MatchPCRight(ch_value,int(ch_label[0])-1,run)
        detector_dictionary[ch_label]=cal_value
    return detector_dictionary

def process_event(line) :
    RootLoad()
    #declare event dictionary
    event_dictionary = dict();
    #split event string
    detector_list = line[1].split(';')
    #loop detectors, put processed object in dictionary
    name_data_pair = detector_list[0].split('?')
    run = -1
    if name_data_pair[0] != 'run' :
        print "Data did not begin with run number!"
        return event_dictionary
    else : run = int(name_data_pair[1])
    for detector in detector_list :
        name_data_pair = detector.split('?')
        if name_data_pair[1] == '' : continue
        elif name_data_pair[0] == 'si-e' : event_dictionary['si-e'] = process_si(name_data_pair[1])
        elif name_data_pair[0] == 'pc-e' : event_dictionary['pc-e'] = process_pc(name_data_pair[1],run)
        elif name_data_pair[0] == 'rf-t' : event_dictionary['rf-t'] = int(name_data_pair[1])
        elif name_data_pair[0] == 'ic-e' : event_dictionary['ic-e'] = int(name_data_pair[1])
    return event_dictionary

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
    for wire in range(1,9) :
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