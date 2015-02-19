#comma sep run list
#run_list  = "048"
run_list  = "009,010,011,012,013"

#file name on hdfs
file_name = "hdfs://gr-gmaster.tamu.edu:9000//data/he8_triumf_0714/tree/carbon_triumf_{{{0}}}*_t.root".format(run_list)
