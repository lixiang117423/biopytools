#!/bin/bash

# 等待 8 小时（8 * 3600 秒 = 28800 秒）
sleep 28800

# 然后执行命令
biopytools kmertools build -i ./02.bam -o 03.kmer_db -t 64
