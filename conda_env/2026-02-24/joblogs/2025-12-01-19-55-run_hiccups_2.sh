mkdir -p hiccups_output_2

java -jar juicer_tools.jar hiccups --cpu \
    -r 1000,2000,5000 \      # 使用更高分辨率：1kb, 2kb, 5kb
    -k KR \
    -f 0.1,0.1,0.1 \
    -p 4,2,1 \               # 峰宽相应调小
    -i 7,5,3 \
    -d 2000,2000,5000 \      # 关键！把距离阈值降下来，允许检测 2kb 以上的 Loop
    aligned/inter_30.hic \
    hiccups_output_2
