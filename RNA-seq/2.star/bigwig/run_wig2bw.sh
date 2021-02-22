
for f in ../*wig
do
    bf=`basename $f`
    tt=`echo $bf | sed -r 's/.wig//g'| sed -r 's/Signal//g' | sed -r 's/.str1//g'| sed -r 's/.out//g'`
    if ! [ -f $tt.bw ]
    then
	echo $tt
	qsub -N $tt -v inp=$f,out=$tt do_wig2bw.sh
    fi
    sleep 1
done


