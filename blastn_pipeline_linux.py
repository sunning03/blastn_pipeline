import os
import argparse
import subprocess
import logging

# 配置日志记录
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
def run_makeblastdb(db_name):
    """
    创建 BLAST 数据库。
    :param db_name: 数据库名称
    """
    # 检查文件是否存在
    fa_file = f"{db_name}.fa"
    fasta_file = f"{db_name}.fasta"
    
    if os.path.exists(fa_file):
        input_file = fa_file
    elif os.path.exists(fasta_file):
        input_file = fasta_file
    else:
        logging.error(f"Neither {fa_file} nor {fasta_file} exists.")
        return
    
    # 构建命令
    command = f"makeblastdb -in {input_file} -dbtype nucl -parse_seqids -out {db_name}"
    logging.info(f"Executing command: {command}")
    
    try:
        result = subprocess.run(command, shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        logging.info("Command output:")
        logging.info(result.stdout)
        if result.stderr:
            logging.warning("Error output:")
            logging.warning(result.stderr)
    except subprocess.CalledProcessError as e:
        logging.error(f"Command execution failed: {e}")
        logging.error(f"Error output: {e.stderr}")
def run_blastn(db_name, query_name):
    """
    运行 BLAST 搜索。
    :param db_name: 数据库名称
    :param query_name: 查询序列文件名
    """

    # 检查文件是否存在
    fa_file = f"{query_name}.fa"
    fasta_file = f"{query_name}.fasta"
    
    if os.path.exists(fa_file):
        query_file = fa_file
    elif os.path.exists(fasta_file):
        query_file = fasta_file
    else:
        logging.error(f"Neither {fa_file} nor {fasta_file} exists.")
        return
    if not os.path.exists(query_file):
        logging.error(f"File {query_file} does not exist.")
        return
    
    output_file = f"{query_name}_{db_name}.blast"
    command = f"blastn -query {query_file} -db {db_name} -evalue 1e-5 -word_size 9 -gapextend 2 -reward 2 -penalty -3 -gapopen 5 -outfmt 6 -num_threads 6 -out {output_file}"
    logging.info(f"Executing command: {command}")
    try:
        result = subprocess.run(command, shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        logging.info("Command output:")
        logging.info(result.stdout)
        if result.stderr:
            logging.warning("Error output:")
            logging.warning(result.stderr)
    except subprocess.CalledProcessError as e:
        logging.error(f"Command execution failed: {e}")
        logging.error(f"Error output: {e.stderr}")

def main():
    """
    主函数，处理命令行参数并运行 BLAST 相关命令。
    """
    parser = argparse.ArgumentParser(description="Run makeblastdb and blastn commands")
    parser.add_argument("db_name", help="Species name of the database")
    parser.add_argument("query_names", nargs="+", help="Names of query sequences")
    
    args = parser.parse_args()
    db_name = args.db_name
    query_names = args.query_names

    run_makeblastdb(db_name)
    for query_name in query_names:
        run_blastn(db_name, query_name)

if __name__ == "__main__":
    main()