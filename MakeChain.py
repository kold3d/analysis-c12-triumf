from ROOT import *
import sys
sys.path.insert(0,"./python")
from runlist import *

#file name on hdfs
chain = TChain("rawData")
list = run_list.split(',')
for run in list :
  chain.Add("/hdfs/data/he8_triumf_0714/tree/he8_triumf_{0}*_t.root".format(run))
