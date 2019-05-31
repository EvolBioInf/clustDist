./clustDist > tmp.out
DIFF=$(diff tmp.out ../data/clustDist.out)
if [ "$DIFF" == "" ] 
then
    printf "Test(clustDist)\tpass\n"
else
    printf "Test(clustDist)\tfail\n"
    echo ${DIFF}
fi

rm tmp.out
