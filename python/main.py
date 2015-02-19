from pyspark import StorageLevel, SparkFiles, SparkContext, SparkConf
import ROOT
import functions as f
from runlist import file_name

def fill_regions(events,scale_factor) :
    s = list()
    for i in range(1,4) : s.append(ROOT.TH1F("s{0}".format(i),"s{0}".format(i),440,1.2,3.4))
    for event in events : 
        value = event.get("cm_energy")
        if value <= 0. : continue
        if value != None :
            label = event.get("si_label")
            if label[0]=='2' : s[0].Fill(value)
            elif (label[0]=='3' and (label[2]=='2' or label[2]=='4')) or \
                 (label[0]=='1' and (label[2]=='1' or label[2]=='3')) : s[1].Fill(value)
            elif (label[0]=='3' and (label[2]=='1' or label[2]=='3')) or \
                 (label[0]=='1' and (label[2]=='2' or label[2]=='4')) : s[2].Fill(value)
    c = ROOT.TCanvas()
    c.Divide(1,3)
    for i in range(1,4) :
        c.cd(i)
        ROOT.Spectra.DivideTargetThickness(s[i-1])
        ROOT.Spectra.CalcSolidAngleNorm(s[i-1],i)
        ##s[i-1].Scale(1./4./3.14159)
        s[i-1].Scale(1./scale_factor)
        s[i-1].Draw()
    out_file = ROOT.TFile("cluster_out.root","recreate")
    c.Write()
    out_file.Close()

if __name__ == "__main__" :
    #scale factor for abs norm
    scale_factor = 1e10
    copyFiles    = True 
    
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
    sc.addFile("Spectra_C.so")
    sc.addFile("cuts.root")
    sc.addFile("dedx_carbon_methane.dat")
    sc.addFile("lookup_table.out")
    sc.addFile("position_cal.txt")
    sc.addFile("wires_scaled_table.out")
    
    #load event file, uses custom hadoop input format
    lines = sc.newAPIHadoopFile(file_name,"edu.tamu.hadoop.RootInputFormat",
                                "org.apache.hadoop.io.IntWritable","org.apache.hadoop.io.Text")

    #fill event dictionaries from root file, process raw events, cut on protons
    proton_events_cm = lines.mapPartitionsWithIndex(lambda x,y : f.process_energy_angle(x,y,copyFiles)) \
                            .flatMap(lambda x : f.get_scattering_events(x)).collect()

    #fill root histogram
    fill_regions(proton_events_cm,scale_factor)
    
