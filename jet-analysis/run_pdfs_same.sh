#!/bin/bash
echo "particle is: " 
echo $1
for i in CT14lo MMHT2014nlo68cl MRST2004qed MRST2004qed_proton CT14nlo MMHT2014lo68cl MRSTMCal MRST2004qed_neutron MRST2007lomod; do
./plot_pdfs_same.py -n $i -q2 $1
done;
