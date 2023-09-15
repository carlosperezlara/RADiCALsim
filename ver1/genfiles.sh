#for X in ele20GeV ele50GeV ele75GeV ele100GeV ele125GeV
for X in ele20GeV ele50GeV ele125GeV pi125GeV
do
    root -b -q -l plotHisto.C\(\"${X}\"\)
    #echo "# 20 GeV electrons (1k)"; cat ${X}.txt > ${X}.txt
done
