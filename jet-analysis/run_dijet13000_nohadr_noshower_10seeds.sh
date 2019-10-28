#!/bin/bash
source /home/petr/.bashrc2
# Energy:
b=13000
#echo "Type number of events:"
c=100000
echo "Running the analysis == rivet-buildplugin RivetMC_DIJET_PB.so MC_DIJET_PB.cc"
rivet-buildplugin RivetMC_DIJET_PB.so MC_DIJET_PB.cc
for s in 61972970 23354883 88615218 70896900 18669684 71917516 25202807 62914023 20426705 23209777;do
echo "Reading file == Herwig read LHC_pb13000.in with seed $s.";
Herwig read LHC_pb13000_no_hadronization.in;
mv LHC.run LHC_${s}.run;
echo "Running Herwig == Herwig run LHC_${s}.run -N ${c} -s ${s}";
Herwig run LHC_${s}.run -N $c -s $s;
done;

