  source ~/.bashrc
  conda activate biopytools

  cd /share/org/YZWL/yzwl_lixg/project/06.longliuxing_BSA/68.三代数据组装和注释/51.两个单倍型的交换图_实验室流程/01.kmer_db

  # 运行导入脚本
  python import_rocksdb.py \
    kmer_matrix.txt.gz \
    kmer_rocksdb \
    --input_delimiter ' ' \
    --batch_size 20000 \
    --bloom_bits 15 \
    --force_overwrite \
    --header_file header.txt \
    --header_db_key rice.32.35

  echo "Import completed at $(date)"
