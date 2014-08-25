from pyspark import SparkFiles, SparkContext, SparkConf
import ROOT
import functions as f

def fill_regions(events,scale_factor) :
    out_file = ROOT.TFile("cluster_out.root","recreate")
    c = ROOT.TCanvas()
    c.Divide(3,2)
    s = list()
    for i in range(1,7) : s.append(ROOT.TH1F("s{0}".format(i),"s{0}".format(i),60,0,3.4))
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
    cuts = list()
    cuts.append(cuts_file.Get("PROTONS_1"))
    cuts.append(cuts_file.Get("PROTONS_2"))
    cuts.append(cuts_file.Get("PROTONS_3"))
    cuts.append(cuts_file.Get("PROTONS_4"))
    cuts.append(cuts_file.Get("PROTONS_5"))
    cuts.append(cuts_file.Get("PROTONS_6"))
    cuts.append(cuts_file.Get("PROTONS_7"))
    cuts.append(cuts_file.Get("PROTONS_8"))
    cuts.append(cuts_file.Get("RF"))
    cuts_file.Close()
    return context.broadcast(cuts)

if __name__ == "__main__" :
    #scale factor for abs norm
    scale_factor = 1.147e9;
    
    #setup cluster
    sconf = SparkConf().setAppName("he8_triumf_analysis")
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
    sc.addFile("dedx_8he_havar.dat")
    sc.addFile("dedx_8he_methane.dat")
    sc.addFile("lookup_table.out")
    sc.addFile("position_cal.txt")
    sc.addFile("wires_scaled_table.out")

    #broadcast cuts to worker nodes
    cuts = export_cuts(sc)

    #comma sep run list
    run_list  = "048,049,052,053,054,055,056,057,058,059,060,061,062,063,064,065,066,067,068,069,070,070,071,073,074,"
    run_list += "075,076,077,078,079,080,081,082,083,084,085,087,088,089,090,091,092,093,094,095,097,098,100,101,"
    run_list += "102,103,104,105,106,107,108"

    file_name = "hdfs://cycdhcp22.tamu.edu:54310/data/he8_triumf_0714/he8_triumf_{{{0}}}001_t.root".format("048")
    print file_name

    #load event file, uses custom hadoop input format
    lines = sc.newAPIHadoopFile(file_name,"edu.tamu.hadoop.RootInputFormat",
                                "org.apache.hadoop.io.IntWritable","org.apache.hadoop.io.Text")

    #fill event dictionaries from root file, process raw events, cut on protons
    proton_events = lines.flatMap(lambda x : f.process_file(x)).map(lambda x : f.process_raw(x)) \
                         .filter(lambda x : f.is_proton(x,cuts.value))

    #lookup cm energy
    proton_events_cm = proton_events.map(lambda x : f.lookup_cm_energy(x)).collect()

    #fill root histogram
    fill_regions(proton_events_cm,scale_factor)
    
