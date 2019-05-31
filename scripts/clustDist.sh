./clustDist ../data/eco29.dist > tmp.out
DIFF=$(diff tmp.out ../data/eco29.nwk)
if [ "$DIFF" == "" ] 
then
    printf "Test(clustDist)\tpass\n"
else
    printf "Test(clustDist)\tfail\n"
    echo ${DIFF}
fi

rm tmp.out
