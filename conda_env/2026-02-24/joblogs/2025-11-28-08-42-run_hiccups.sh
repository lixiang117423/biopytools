# 确保输出目录存在（虽然 Juicer 通常会自动创建，但为了保险）
mkdir -p hiccups_output

# 运行命令
java -Xmx128g -jar ~/software/juicer/scripts/juicer_tools.jar hiccups --cpu \
    -r 5000,10000 \
    -k KR \
    -f 0.1,0.1 \
    -p 4,2 \
    -i 7,5 \
    -d 20000,20000 \
    aligned/inter_30.hic \
    hiccups_output