#!/bin/bash
echo "Type Energy:"
read b
#echo "Type LEP.yoda for herwig and Rivet.yoda for python"
a="LHC.yoda"
f="Rivet.yoda"
echo "Type number of events:"
read c
echo "Running the analysis == rivet-buildplugin RivetMC_DIJET_PB.so MC_DIJET_PB.cc"
rivet-buildplugin RivetMC_DIJET_PB.so MC_DIJET_PB.cc
if [ $b == "13000" ];
then
echo "Reading file == Herwig read LHC_pb13000.in"
Herwig read LHC_pb13000.in
echo "Running Herwig == Herwig run LHC.run -N ${c}"
Herwig run LHC.run -N $c;
fi
if [ $b == "900" ];
then
echo "Reading file == Herwig read LHC_pb900.in"
Herwig read LHC_pb900.in
echo "Running Herwig == Herwig run LHC.run -N ${c}"
Herwig run LHC.run -N $c;
fi


