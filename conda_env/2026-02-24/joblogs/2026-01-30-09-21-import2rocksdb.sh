  # 运行导入脚本
  python3 import_rocksdb.py \
    kmer_matrix.txt.gz \
    kmer_rocksdb \
    --input_delimiter ' ' \
    --batch_size 20000 \
    --bloom_bits 15 \
    --force_overwrite \
    --header_file header.txt \
    --header_db_key rice.32.35

  echo "Import completed at $(date)"
