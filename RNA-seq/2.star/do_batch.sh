

fqdir='fqs.p33'

# PE reads first, as they are likely to be faster.

for f in  ../$fqdir/*.fq.gz
do
    root=`basename $f`
    
    if [[ $root == *.p1.fq* ]]
    then
	echo $root | grep -q "Mmmm"
	if [ $? -eq 0 ]; then
		echo $root
		echo "pass"
	else
        	bf=`echo $root | sed -r 's/.p1.fq//g' | sed -r 's/.gz//g'`
        	p1=`echo $root`
        	p2=`echo $root | sed 's/\.p1/.p2/g'`
        	tt=`echo $bf.results.txt`
        	if ! [ -f $tt ]
        	then
	    		#echo $p1,$p2
            		echo PE ... $tt
            		qsub -N $bf.star.pe -v p1=../$fqdir/$p1,p2=../$fqdir/$p2,out=$bf do_pe.sh
            		sleep 2
        	fi
	fi
    else # not a PE read
        if [[ $root != *.p2.fq* ]]
        then
            bf=`echo $root | sed -r 's/.fq//g' | sed -r s/.gz//g`
            tt=`echo $bf.results.txt`
            if ! [ -f $tt ]
            then
                echo SE ... $tt
                qsub -N $bf.star.se -v inp=../$fqdir/$root,out=$bf do_se.sh
                sleep 2
            fi
        fi
    fi
done


