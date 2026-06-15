import os
from tqdm import tqdm
import pandas as pd


class Record(object):
    '''
    One line information in vcf file
    '''

    def __init__(self, line):
        info = line.split("\t")
        self.line = line
        self.CHROM = info[0]
        self.POS = info[1]
        self.ID = info[2]
        self.REF = info[3]
        self.ALT = info[4]
        self.QUAL = info[5]
        self.FILTER = info[6]
        self.INFO = [{pair_lst[0]: pair_lst[1] if len(pair_lst) > 1 else ""} for pair_lst in
                     [pair.split("=") for pair in info[7].split(";")]]
        self.FORMAT = info[8].split(":")
        self.sample_num = len(info) - 9
        self.GT = []
        for i in range(self.sample_num):
            GT_value = info[8 + i + 1].split(":")
            GT_dict = {}
            for g in range(len(GT_value)):
                GT_dict[self.FORMAT[g]] = GT_value[g]
            self.GT.append(GT_dict)


class VCF(object):
    '''
    VCF class, read VCF, write VCF, get VCF information
    '''

    def __init__(self, uncompress_vcf):
        self.header = []
        self.reader = open(uncompress_vcf, 'r')
        self.line = self.reader.readline().strip()
        while self.line.startswith('#'):
            self.header.append(self.line)
            self.line = self.reader.readline().strip()
        self.record = Record(self.line)

    def __iter__(self):
        return self

    def __next__(self):
        self.line = self.reader.readline().strip()
        if self.line != "":
            self.record = Record(self.line)
            return self.record
        else:
            self.reader.close()
            raise StopIteration()

    def reader_close(self):
        self.reader.close()


class VCF2Excel(object):

    def __init__(self, file_path, file_name, save_dir):
        super(VCF2Excel, self).__init__()
        self.file_path = file_path
        self.file_name = file_name
        self.save_dir = save_dir

    # def run(self):
    #     vcf = VCF(self.file_path)
    #     data = []
    #     for r in tqdm(vcf, desc="vcf to excel"):
    #         try:
    #             row = []
    #             row.append(str(r.CHROM))
    #             row.append(int(r.POS))
    #             row.append(r.REF)
    #             row.append(r.ALT)
    #             for i in r.GT:
    #                 row.append(int(i["AD"].split(",")[0]))
    #                 row.append(int(i["AD"].split(",")[1]))
    #             data.append(row)
    #         except:
    #             pass

    #     df = pd.DataFrame(data)
    #     save_path = os.path.join(self.save_dir, self.file_name + ".csv")
    #     df.to_csv(save_path, header=False, index=False)

    #     return save_path

    def run(self):
        vcf = VCF(self.file_path)
        data = []
        error_count = 0
        
        for r in tqdm(vcf, desc="vcf to excel"):
            try:
                row = []
                row.append(str(r.CHROM))
                row.append(int(r.POS))
                row.append(r.REF)
                row.append(r.ALT)
                for i in r.GT:
                    ad_values = i["AD"].split(",")
                    if len(ad_values) >= 2:  # 确保有足够的数据
                        row.append(int(ad_values[0]))
                        row.append(int(ad_values[1]))
                    else:
                        raise ValueError("AD field incomplete")
                
                # 只添加完整的行
                if len(row) >= 4:  # 至少有 CHROM, POS, REF, ALT
                    data.append(row)
            except Exception as e:
                error_count += 1
                if error_count <= 5:  # 只打印前5个错误
                    print(f"Warning: Skipping malformed line - {e}")
                continue
        
        print(f"Total lines processed: {len(data)}, Errors: {error_count}")
        
        df = pd.DataFrame(data)
        save_path = os.path.join(self.save_dir, self.file_name + ".csv")

        # 明确指定CSV格式（移除lineterminator以兼容旧版pandas）
        # Specify CSV format (remove lineterminator for compatibility with older pandas)
        df.to_csv(save_path, header=False, index=False, sep=',')

        print(f"CSV saved to: {save_path}")
        print(f"CSV shape: {df.shape}")

        return save_path