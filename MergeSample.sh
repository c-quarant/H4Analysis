#!bin/bash

output="btf2016_Rm_8m"
outputMerge="$output"
outputMerge+="_merged_2460_2464.root"
dir="eos/cms/store/group/dpg_ecal/comm_ecal/upgrade/testbeam/TimingTB_T9_Summer2016/TimingTB_BTF_June2016/ntuples/v1"
hadd="python macros/merger.py $dir/$outputMerge"

for run in 2460 2461 2462 2463 2464; do 
    outputRun="$output"
    outputRun+="_$run"
    hadd+=" $dir/$outputRun.root"; 
done

$hadd
