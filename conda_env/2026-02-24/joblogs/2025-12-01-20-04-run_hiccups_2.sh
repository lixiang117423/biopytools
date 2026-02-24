mkdir -p hiccups_output_2

java -Xmx128g -jar ~/software/juicer/scripts/juicer_tools.3.0.0.jar hiccups --cpu \
    -r 1000,2000,5000 -k KR -f 0.1,0.1,0.1 -p 4,2,1 -i 7,5,3 -d 2000,2000,5000 aligned/inter_30.hic hiccups_output_2
