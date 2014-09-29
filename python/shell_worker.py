import ROOT

def select_hist(command,conditions,filename) :
    hist = None
    ROOT.gROOT.SetBatch()
    file = ROOT.THDFSFile(filename[1])
    tree = file.Get("rawData")
    events = tree.Draw(command,conditions)
    if events > 0 : 
        hist = ROOT.gPad.GetPrimitive("h").Clone(command)
        hist.SetDirectory(0)
    file.Close()
    return hist
    
def select_hist_1d(command,conditions,binsx,lowx,highx,filename) :
    new_command = "{0}>>h({1},{2},{3})".format(command,binsx,lowx,highx)
    return select_hist(new_command,conditions,filename)

def select_hist_2d(command,conditions,binsx,lowx,highx,binsy,lowy,highy,filename) :
    new_command = "{0}>>h({1},{2},{3},{4},{5},{6})".format(command,binsy,lowy,highy,binsx,lowx,highx)
    return select_hist(new_command,conditions,filename)
