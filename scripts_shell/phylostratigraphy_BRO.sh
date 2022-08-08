for i in $(cat $1); do 

og=$(grep $i'.*_BRO' Orthogroup_Sequences/*fa | tail -1 | awk -F ":" '{print $1}');

echo -n $og $i" " >> $1.count.tmp

if grep -q PHY $og; then 

if cat $og | tr -d "\n" | grep -qe "BGM.*BAT" -qe "BAT.*BGM"; then echo ubiquitous >> $1.count.tmp; 

else echo shared_with_phyllium >> $1.count.tmp;

fi;

elif cat $og | tr -d "\n" | grep -qe "BGM.*BAT" -qe "BAT.*BGM"; then echo ubiquitous_in_bacillus >> $1.count.tmp;

elif grep -qE '(BGM|BAT)' $og; then echo shared_in_bacillus >> $1.count.tmp; 

else echo species_specific >> $1.count.tmp;

fi

unset og

done

ubiquitous=$(grep -wc ubiquitous $1.count.tmp)
shared_with_phyllium=$(grep -wc "shared_with_phyllium" $1.count.tmp)
ubiquitous_in_bacillus=$(grep -wc "ubiquitous_in_bacillus" $1.count.tmp)
shared_in_bacillus=$(grep -wc "shared_in_bacillus" $1.count.tmp)
species_specific=$(grep -wc "species_specific" $1.count.tmp)

echo -e "total genes \t \t \t $(wc -l $1.count.tmp)" >> $1.result
echo -e "ubiquitous \t \t $ubiquitous" >> $1.result
echo -e "shared with phyllium \t $shared_with_phyllium" >> $1.result
echo -e "ubiquitous in bacillus \t $ubiquitous_in_bacillus" >> $1.result
echo -e "shared in bacillus \t $shared_in_bacillus" >> $1.result
echo -e "species specific \t $species_specific" >> $1.result

grep specific $1.count.tmp > $1.species_specific_genes.result

mv $1.count.tmp $1.count