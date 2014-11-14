import ROOT
from runlist import file_name
import shell_worker as sw

def add_hists_1d(hists,command,binsx,lowx,highx) :
    sum_hist = ROOT.TH1F(command,command,binsx,lowx,highx)
    for i in xrange(0,len(hists)) :
        if hists[i] != None : sum_hist.Add(hists[i],1.)
    return sum_hist

def draw_1d(sc,command,conditions,binsx,lowx,highx) :
    lines = sc.newAPIHadoopFile(file_name,"edu.tamu.hadoop.RootInputFormat",
                                "org.apache.hadoop.io.IntWritable","org.apache.hadoop.io.Text")
    hists = lines.map(lambda x : sw.select_hist_1d(command,conditions,binsx,lowx,highx,x)).collect()
    sum_hist =  add_hists_1d(hists,command,binsx,lowx,highx)
    sum_hist.Draw()
    return sum_hist

def add_hists_2d(hists,command,binsx,lowx,highx,binsy,lowy,highy) :
    sum_hist = ROOT.TH2F(command,command,binsx,lowx,highx,binsy,lowy,highy)
    for i in xrange(0,len(hists)) :
        if hists[i] != None : sum_hist.Add(hists[i],1.)
    return sum_hist

def draw_2d(sc,command,conditions,binsx,lowx,highx,binsy,lowy,highy) :
    lines = sc.newAPIHadoopFile(file_name,"edu.tamu.hadoop.RootInputFormat",
                                "org.apache.hadoop.io.IntWritable","org.apache.hadoop.io.Text")
    hists = lines.map(lambda x : sw.select_hist_2d(command,conditions,binsx,lowx,highx,binsy,lowy,highy,x)).collect()
    sum_hist =  add_hists_2d(hists,command,binsx,lowx,highx,binsy,lowy,highy)
    sum_hist.Draw()
    return sum_hist
