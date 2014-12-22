#comma sep run list
#run_list  = "048"
run_list  = "048,049,052,053,054,055,056,057,063,064,065,066,067,068,069,070,071,073,074,"
run_list += "075,076,077,079,080,081,082,083,084,085,087,088,089,090,091,092,093,094,095,096,097,098,100,101,"
run_list += "102,103,104,105,106,107"

#file name on hdfs
file_name = "hdfs://gr-gmaster.tamu.edu:9000//data/he8_triumf_0714/tree/he8_triumf_{{{0}}}*_t.root".format(run_list)
