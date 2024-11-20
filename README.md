首先需要安装blast，可参考教程[blast | Sunning's Blog (sunning03.github.io)](https://sunning03.github.io/2024/03/12/blast/)

## Linux版使用



`python blastn_pipeline_linux.py -h`

```
positional arguments:
  db_name      Species name of the database
  query_names  Names of query sequences

optional arguments:
  -h, --help   show this help message and exit
```

`python blastn_pipeline_linux.py db_name query_name_1 query_name_2 ……query_name_n`  

1. db_name为makeblastdb的物种名
2. query_names为需要比对的物种名（可一次输入多个）

使用案例：

`python blast_pipeline_linux.py Pru_kan Pru_yed Pru_sch`

其中，Pru_kan为建库物种，Pru_yed和Pru_sch是需要比对的物种

三个物种的fa文件已放在Example文件夹内。

## Note

1. 目前该脚本仅支持后缀为.fa的文件，请修改自己的文件后缀以适用脚本

2. makeblastdb和blastn的参数设置如下，如果需要修改参数，请自行修改脚本

   ```
   makeblastdb -in {fa_file} -dbtype nucl -parse_seqids -out {db_name}
   blastn -query {query_file} -db {db_name} -evalue 1e-5 -word_size 9 -gapextend 2 -reward 2 -penalty -3 -gapopen 5 -outfmt 6 -num_threads 6 -out {output_file}
   ```

   