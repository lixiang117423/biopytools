# 自动识别多种压缩格式
import rocksdb
import gzip
import os
import shutil
import argparse
import traceback

def get_rocksdb_options(bloom_bits_per_key=15):
    """
    Configures RocksDB options for efficient bulk loading.
    NOTE: This version no longer handles compression fallback. The calling function does.
    """
    opts = rocksdb.Options()
    opts.create_if_missing = True
    
    opts.table_factory = rocksdb.BlockBasedTableFactory(
        filter_policy=rocksdb.BloomFilterPolicy(bits_per_key=bloom_bits_per_key),
        block_cache=rocksdb.LRUCache(128 * (1024**2)),
    )
    
    cpu_cores = os.cpu_count()
    if cpu_cores and cpu_cores > 1:
        if cpu_cores <= 2:
            opts.max_background_compactions = 1
            opts.max_background_flushes = 1
        else:
            num_total_bg_threads = max(2, cpu_cores // 2)
            if num_total_bg_threads < 4 and cpu_cores > 2:
                num_total_bg_threads = min(4, cpu_cores)
            opts.max_background_flushes = max(1, num_total_bg_threads // 4 if num_total_bg_threads > 3 else 1)
            opts.max_background_compactions = max(1, num_total_bg_threads - opts.max_background_flushes)
    else:
        opts.max_background_compactions = 1
        opts.max_background_flushes = 1
        
    opts.write_buffer_size = 128 * (1024**2)
    opts.max_write_buffer_number = 4
    opts.target_file_size_base = 128 * (1024**2)
    opts.max_bytes_for_level_base = 512 * (1024**2)
    opts.level0_file_num_compaction_trigger = 10
    opts.level0_slowdown_writes_trigger = 20
    opts.level0_stop_writes_trigger = 30
    
    return opts


def import_gz_to_rocksdb(gz_file_path, db_path, input_delimiter='\t',
                         value_join_char='',
                         batch_size=10000, bloom_bits=15, force_overwrite=False,
                         header_file_path=None, header_db_key="kmer_header",
                         header_storage_delimiter='\t'):
    print(f"Starting import from '{gz_file_path}' to RocksDB at '{db_path}'...")
    if header_file_path:
        print(f"INFO: Header file: '{header_file_path}', DB key: '{header_db_key}', Storage delimiter for header fields: '{header_storage_delimiter}'")
    print(f"Input delimiter: '{input_delimiter}'")
    print(f"Data value join character (for storage): '{value_join_char}' (Empty means concatenation)")

    if os.path.exists(db_path):
        if force_overwrite:
            print(f"WARNING: Overwriting existing DB at '{db_path}'.")
            shutil.rmtree(db_path)
        else:
            print(f"ERROR: DB path '{db_path}' exists. Use --force_overwrite.")
            return False

    db = None
    compression_preferences = [
        ("ZSTD", rocksdb.CompressionType.zstd_compression),
        ("Snappy", rocksdb.CompressionType.snappy_compression),
        ("LZ4", rocksdb.CompressionType.lz4_compression),
        ("None", rocksdb.CompressionType.no_compression)
    ]

    for comp_name, comp_type in compression_preferences:
        print(f"INFO: Attempting to open DB with {comp_name} compression...")
        rocks_opts = get_rocksdb_options(bloom_bits_per_key=bloom_bits)
        
        try:
            # AttributeError can happen if the compression type constant doesn't exist at all
            rocks_opts.compression = comp_type
            db = rocksdb.DB(db_path, rocks_opts)
            print(f"SUCCESS: Opened DB with {comp_name} compression.")
            break 
        # --- 这是关键的修改 ---
        except (rocksdb.errors.InvalidArgument, AttributeError) as e:
        # ----------------------
            if "not linked with the binary" in str(e) or "has no attribute" in str(e):
                print(f"INFO: {comp_name} is not supported. Trying next option.")
                if os.path.exists(db_path):
                    shutil.rmtree(db_path)
            else:
                print(f"Error opening/creating RocksDB with {comp_name}: {e}")
                return False
    
    if not db:
        print("ERROR: Failed to open RocksDB with any available compression type.")
        return False

    batch = rocksdb.WriteBatch()
    count_in_batch = 0
    imported_records_count = 0
    physical_lines_read = 0
    current_line_idx = -1

    header_fields_from_file = []
    if header_file_path:
        if os.path.exists(header_file_path):
            try:
                with open(header_file_path, 'rt', encoding='utf-8') as hf:
                    for line in hf:
                        field_name = line.strip()
                        if field_name:
                            header_fields_from_file.append(field_name)
                
                if header_fields_from_file:
                    header_to_store_in_db = header_storage_delimiter.join(header_fields_from_file)
                    key_bytes = header_db_key.encode('utf-8')
                    value_bytes = header_to_store_in_db.encode('utf-8')
                    batch.put(key_bytes, value_bytes)
                    count_in_batch += 1
                    print(f"INFO: Added processed header to batch (key='{header_db_key}'): '{header_to_store_in_db}'")
                else:
                    print(f"WARNING: Header file '{header_file_path}' was empty.")
            except Exception as e:
                print(f"ERROR processing header file '{header_file_path}': {e}")
        else:
            print(f"WARNING: Header file '{header_file_path}' not found.")

    try:
        open_func = gzip.open if gz_file_path.endswith('.gz') else open
        
        with open_func(gz_file_path, 'rt', encoding='utf-8') as infile:
            for current_line_idx, line in enumerate(infile, 0):
                physical_lines_read = current_line_idx + 1
                line_content = line.rstrip('\n')
                if not line_content: continue

                parts = line_content.split(input_delimiter, 1)
                key_str = parts[0]
                
                if not key_str: continue

                value_str_combined = ""
                if len(parts) > 1:
                    remaining_value_part = parts[1]
                    value_columns = remaining_value_part.split(input_delimiter)
                    value_str_combined = value_join_char.join(value_columns)
                                
                key_bytes = key_str.encode('utf-8')
                value_bytes = value_str_combined.encode('utf-8')
                batch.put(key_bytes, value_bytes)
                count_in_batch += 1

                if count_in_batch >= batch_size:
                    db.write(batch)
                    imported_records_count += batch.count()
                    if physical_lines_read % (batch_size * 100) == 0:
                        print(f"Processed lines: {physical_lines_read}, DB records imported so far: {imported_records_count}")
                    batch.clear()
                    count_in_batch = 0
            
            if batch.count() > 0:
                db.write(batch)
                imported_records_count += batch.count()

        print("-" * 30)
        print("Import finished successfully!")
        data_records_count = imported_records_count - (1 if header_fields_from_file else 0)
        print(f"Total lines read from input file: {physical_lines_read}")
        print(f"Total data records imported: {data_records_count}")
        if header_fields_from_file:
            print("Header record was also imported.")
        return True

    except Exception as e:
        err_line_msg = f"at input file line {current_line_idx + 1}" if current_line_idx != -1 else "during processing"
        print(f"ERROR: An unexpected error occurred {err_line_msg}: {e}")
        traceback.print_exc()
    finally:
        pass
    return False

def main():
    parser = argparse.ArgumentParser(description="Import data from a file into RocksDB.")
    parser.add_argument("input_file", help="Input data file (can be .gz or plain text).")
    parser.add_argument("rocksdb_path", help="RocksDB output path.")
    parser.add_argument("--input_delimiter", default='\t', help="Delimiter in the input file (default: tab). Use ' ' for space.")
    parser.add_argument("--batch_size", type=int, default=20000, help="Number of records to write in a single batch.")
    parser.add_argument("--bloom_bits", type=int, default=15, help="Bits per key for the Bloom filter.")
    parser.add_argument("--force_overwrite", action="store_true", help="Overwrite the existing RocksDB database.")
    parser.add_argument("--header_file", help="Path to a header file with one column name per line.")
    parser.add_argument("--header_db_key", default="kmer_header", help="The database key for the header string.")
    
    args = parser.parse_args()

    if not os.path.exists(args.input_file):
        print(f"ERROR: Input file '{args.input_file}' does not exist.")
        return

    success = import_gz_to_rocksdb(
        args.input_file,
        args.rocksdb_path,
        input_delimiter=args.input_delimiter,
        batch_size=args.batch_size,
        bloom_bits=args.bloom_bits,
        force_overwrite=args.force_overwrite,
        header_file_path=args.header_file,
        header_db_key=args.header_db_key,
    )
    if success:
        print("\nImport completed successfully.")
    else:
        print("\nImport failed.")

if __name__ == "__main__":
    main()
