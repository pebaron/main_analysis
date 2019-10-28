#!/bin/bash
source /home/petr/.bashrc2
# Energy:
b=900
#echo "Type number of events:"
c=100000
echo "Running the analysis == rivet-buildplugin RivetMC_DIJET_PB.so MC_DIJET_PB.cc"
rivet-buildplugin RivetMC_DIJET_PB.so MC_DIJET_PB.cc
for s in 61972970 23354883 88615218 70896900 18669684 71917516 25202807 62914023 20426705 23209777;do
echo "Reading file == Herwig read LHC_pb900.in with seed $s.";
Herwig read LHC_pb900.in;
mkdir ${s}_900
mv LHC.run ./${s}_900/LHC_${s}_900.run;
echo "Running Herwig == Herwig run LHC_${s}_900.run -N ${c} -s ${s}";
Herwig run ./${s}_900/LHC_${s}_900.run -N $c -s $s;
done;

