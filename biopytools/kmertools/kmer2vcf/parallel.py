"""
多线程处理模块|Multithreading Processing Module
"""

import os
import queue
import threading
from typing import Callable, Any, Tuple
from .utils import format_number, open_input


class ChunkProcessor:
    """分块处理器|Chunk Processor for streaming multithreading"""

    def __init__(self, config, logger):
        """
        初始化分块处理器|Initialize chunk processor

        Args:
            config: 配置对象|Configuration object
            logger: 日志器|Logger object
        """
        self.config = config
        self.logger = logger
        self.chunk_size = config.chunk_size
        self.num_threads = config.threads

    def process_file_in_chunks(
        self,
        input_file: str,
        output_file: str,
        process_func: Callable,
        header_line: str = None,
        skip_header: bool = True
    ) -> int:
        """
        分块读取文件并多线程处理|Read file in chunks and process with multithreading

        Args:
            input_file: 输入文件|Input file path
            output_file: 输出文件|Output file path
            process_func: 处理函数|Processing function (takes list of lines, returns processed result)
            header_line: 表头行|Header line to write
            skip_header: 是否跳过输入文件表头|Whether to skip input file header

        Returns:
            int: 总行数|Total lines processed
        """
        self.logger.info(f"多线程处理文件|Multithreaded processing: {self.num_threads} threads")
        self.logger.info(f"块大小|Chunk size: {format_number(self.chunk_size)} lines")

        # 创建任务队列和结果队列|Create task queue and result queue
        # 限制队列深度以控制内存|Limit queue depth to control memory
        max_queue_size = 5
        task_queue = queue.Queue(maxsize=max_queue_size)
        result_queue = queue.Queue(maxsize=max_queue_size)

        # 启动worker线程|Start worker threads
        workers = []
        for i in range(self.num_threads):
            worker = threading.Thread(
                target=self._worker,
                args=(task_queue, result_queue, process_func, i),
                daemon=True
            )
            worker.start()
            workers.append(worker)

        # 使用共享变量保存producer结果|Use shared variable for producer results
        producer_result = {'total_lines': 0, 'header_line': header_line, 'chunk_id': 0}

        def producer():
            nonlocal skip_header
            total_lines = 0
            chunk_id = 0
            local_header = header_line

            with open_input(input_file) as f:
                if skip_header and local_header is None:
                    # 读取并保存header|Read and save header
                    first_line = f.readline()
                    if first_line.startswith('#') or '\t' in first_line or ' ' in first_line:
                        local_header = first_line.strip()
                    else:
                        # 第一行不是header，需要放回|First line is not header, put it back
                        f.seek(0)
                        skip_header = False

                chunk = []
                for line in f:
                    total_lines += 1
                    chunk.append(line)

                    if len(chunk) >= self.chunk_size:
                        task_queue.put((chunk_id, chunk))
                        chunk_id += 1
                        chunk = []

                        # 进度显示|Progress display
                        if total_lines % (self.chunk_size * 10) == 0:
                            self.logger.info(f"已读取|Read: {format_number(total_lines)} lines")

                # 处理剩余的行|Process remaining lines
                if chunk:
                    task_queue.put((chunk_id, chunk))
                    chunk_id += 1

            # 发送结束信号|Send end signal
            for _ in range(self.num_threads):
                task_queue.put(None)

            # 保存结果到共享变量|Save results to shared variable
            producer_result['total_lines'] = total_lines
            producer_result['header_line'] = local_header
            producer_result['chunk_id'] = chunk_id

        # 启动生产者|Start producer
        producer_thread = threading.Thread(target=producer, daemon=True)
        producer_thread.start()

        # 等待producer完成以获取总chunk数|Wait for producer to get total chunks
        producer_thread.join()
        total_chunks = producer_result['chunk_id']

        # 消费者：收集结果并按顺序写入|Consumer: collect results and write in order
        total_lines = self._write_results_sequentially(
            result_queue, output_file, producer_result['header_line'], total_chunks
        )

        # 等待所有worker完成|Wait for all workers to complete
        for worker in workers:
            worker.join()

        return total_lines

    def _worker(
        self,
        task_queue: queue.Queue,
        result_queue: queue.Queue,
        process_func: Callable,
        worker_id: int
    ):
        """
        Worker线程处理函数|Worker thread processing function

        Args:
            task_queue: 任务队列|Task queue
            result_queue: 结果队列|Result queue
            process_func: 处理函数|Processing function
            worker_id: Worker ID
        """
        while True:
            task = task_queue.get()

            if task is None:
                # 结束信号|End signal
                break

            chunk_id, chunk = task

            try:
                # 处理数据块|Process chunk
                result = process_func(chunk)
                result_queue.put((chunk_id, result))
            except Exception as e:
                self.logger.error(f"Worker {worker_id} 处理块 {chunk_id} 时出错|Error processing chunk {chunk_id}: {e}")
                result_queue.put((chunk_id, None))

            task_queue.task_done()

    def _write_results_sequentially(
        self,
        result_queue: queue.Queue,
        output_file: str,
        header_line: str,
        total_chunks: int
    ) -> int:
        """
        按顺序写入结果|Write results in sequential order

        Args:
            result_queue: 结果队列|Result queue
            output_file: 输出文件|Output file
            header_line: 表头|Header line
            total_chunks: 总块数|Total number of chunks

        Returns:
            int: 总行数|Total lines written
        """
        # 使用缓冲区来保持顺序|Use buffer to maintain order
        pending_chunks = {}
        next_chunk_id = 0
        total_lines = 0

        with open(output_file, 'w') as f_out:
            # 写入表头|Write header
            if header_line:
                f_out.write(header_line + '\n')

            # 收集并按顺序写入|Collect and write in order
            received_chunks = 0

            while received_chunks < total_chunks:
                chunk_id, result = result_queue.get()

                if result is None:
                    # 处理失败，跳过|Processing failed, skip
                    received_chunks += 1
                    continue

                received_chunks += 1

                # 如果是期望的chunk，写入；否则缓存|If expected chunk, write; otherwise buffer
                if chunk_id == next_chunk_id:
                    f_out.write(result)
                    total_lines += result.count('\n')

                    # 检查并写入缓存的chunks|Check and write buffered chunks
                    next_chunk_id += 1
                    while next_chunk_id in pending_chunks:
                        f_out.write(pending_chunks[next_chunk_id])
                        total_lines += pending_chunks[next_chunk_id].count('\n')
                        del pending_chunks[next_chunk_id]
                        next_chunk_id += 1
                else:
                    # 缓存乱序到达的chunk|Buffer out-of-order chunk
                    pending_chunks[chunk_id] = result

        return total_lines


def convert_chunk_to_binary(chunk: list) -> str:
    """
    转换数据块为0/1格式（Worker处理函数）|Convert chunk to binary format (worker processing function)

    Args:
        chunk: 行列表|List of lines

    Returns:
        str: 转换后的字符串|Converted string
    """
    result_lines = []

    for line in chunk:
        line = line.rstrip('\n')

        # 尝试制表符分隔|Try tab delimiter first
        parts = line.split('\t')

        # 如果制表符分隔结果少于2列，尝试空格分隔|If tab gives <2 columns, try space
        if len(parts) < 2:
            parts = line.split()

        if len(parts) < 2:
            continue

        kmer_id = parts[0]
        abundances = parts[1:]

        # 验证丰度数据是否为数字（跳过header行）|Validate abundances are numbers (skip header lines)
        try:
            # 尝试转换第一个值为整数来验证|Try converting first value to validate
            int(abundances[0])
        except (ValueError, IndexError):
            # 如果不是数字，跳过这一行（可能是header）|If not a number, skip this line (likely header)
            continue

        # 转换丰度为0/1|Convert abundances to binary
        binary_values = ['1' if int(abd) > 0 else '0' for abd in abundances]

        # 写入临时文件|Write to temporary file
        result_lines.append(f"{kmer_id}\t{' '.join(binary_values)}\n")

    return ''.join(result_lines)
