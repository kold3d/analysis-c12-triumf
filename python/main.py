from pyspark import SparkFiles, SparkContext, SparkConf
import ROOT
import functions as f

def fill_regions(events,scale_factor) :
    out_file = ROOT.TFile("cluster_out.root","recreate")
    c = ROOT.TCanvas()
    c.Divide(3,2)
    s = list()
    for i in range(1,7) : s.append(ROOT.TH1F("s{0}".format(i),"s{0}".format(i),100,0,3.))
    for event in events : 
        value = event.get("cm_energy")
        if value <= 0. : continue
        if value != None :
            which = 0
            if event.get("pc_wire") > 5 : which = 1
            bounds = list()
            for i in range(0,6) : bounds.append(ROOT.Spectra.LookupPCBound(which,i+1,value))
            position = abs(event.get("pc_position"))
            label = event.get("si_label")
            if label[0]=='2' and position <= bounds[0].second : s[0].Fill(value)
            elif label[0]=='2' and position > bounds[1].first : s[1].Fill(value)
            elif (label[0]=='3' or label[0]=='1') and position < bounds[2].second : s[2].Fill(value)
            elif (label[0]=='3' or label[0]=='1') and position > bounds[3].first and \
                 position <= bounds[3].second : s[3].Fill(value)
            elif (label[0]=='3' or label[0]=='1') and position > bounds[4].first and \
                 position <= bounds[4].second : s[4].Fill(value)
            elif (label[0]=='3' or label[0]=='1') and position > bounds[5].first : s[5].Fill(value)
    for i in range(1,7) :
        c.cd(i)
        ROOT.Spectra.DivideTargetThickness(s[i-1])
        ROOT.Spectra.CalcSolidAngleFast(s[i-1],i)
        s[i-1].Scale(1./scale_factor)
        s[i-1].Draw()
    c.Write()
    out_file.Close()

#function opens cut files, reads them cuts, and broadcasts cuts to worker nodes
def export_cuts(context) : 
    cuts_file = ROOT.TFile("cuts.root","read")
    p_cuts = list()
    p_cuts.append(cuts_file.Get("PROTONS_1"))
    p_cuts.append(cuts_file.Get("PROTONS_2"))
    p_cuts.append(cuts_file.Get("PROTONS_3"))
    p_cuts.append(cuts_file.Get("PROTONS_4"))
    p_cuts.append(cuts_file.Get("PROTONS_5"))
    p_cuts.append(cuts_file.Get("PROTONS_6"))
    p_cuts.append(cuts_file.Get("PROTONS_7"))
    p_cuts.append(cuts_file.Get("PROTONS_8"))
    cuts_file.Close()
    return context.broadcast(p_cuts)

if __name__ == "__main__" :
    #scale factor for abs norm
    scale_factor = 1.e10;
    
    #setup cluster
    sconf = SparkConf().setAppName("c1_triumf_analysis")
    sc = SparkContext(conf=sconf)
   
    #build root libraries
    ROOT.gROOT.SetBatch()
    ROOT.gSystem.CompileMacro("EnergyLoss.C")
    ROOT.gSystem.CompileMacro("Calibrations.C")
    ROOT.Calibrations.InitParameters()
    ROOT.gSystem.CompileMacro("Spectra.C")
    ROOT.Spectra.ReadPCBoundTable()
    ROOT.Spectra.ReadSolidAngleTable()
    ROOT.gSystem.CompileMacro("EnergyAngle.C")

    #Add Analysis Files
    sc.addFile("EnergyLoss_C.so")
    sc.addFile("Calibrations_C.so")
    sc.addFile("EnergyAngle_C.so")
    sc.addFile("dEdx_carbon_methane_290K_400torr.dat")
    sc.addFile("dEdx_proton_methane_290K_400torr.dat")
    sc.addFile("lookup_table.out")
    sc.addFile("position_cal.txt")

    #broadcast cuts to worker nodes
    p_cuts = export_cuts(sc)

    #load event file
    lines = sc.textFile("hdfs://cycdhcp22.tamu.edu:54310/data/c12_triumf_0714/carbon_triumf_*_ascii.txt")

    #fill event dictionaries from ascii, find cm energies
    events = lines.map(lambda x : f.process_event(x))
    ea_events = events.map(lambda x : f.process_energy_angle(x))
    proton_ea_events = ea_events.filter(lambda x : f.is_proton(x,p_cuts.value)).collect()

    #fill root histogram
    fill_regions(proton_ea_events,scale_factor)
    
