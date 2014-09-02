import ROOT

def select_hist_1d(command,conditions,binsx,lowx,highx,filename) :
    hist = None
    ROOT.gROOT.SetBatch()
    file = ROOT.THDFSFile(filename[1])
    tree = file.Get("rawData")
    events = tree.Draw("{0}>>h({1},{2},{3})".format(command,binsx,lowx,highx),conditions)
    if events > 0 : 
        hist = ROOT.gPad.GetPrimitive("h").Clone(command)
        hist.SetDirectory(0)
    file.Close()
    return hist
