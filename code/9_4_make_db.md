# 使用 eggnog-mapper 创建 GO 注释数据库

在 http://eggnog-mapper.embl.de/ 上传 Roary 输出的 `data/roary_result/pan_genome_reference.fa` 文件；数据类型选择 `CDS`，**Annotation options** 选择 `Brucellaceae` (布鲁氏菌科)，否则注释结果为空；

在线运行给出的指令是 

```shell
emapper.py --cpu 20 --mp_start_method forkserver --data_dir /dev/shm/ -o out --output_dir /emapper_web_jobs/emapper_jobs/user_data/MM_wpfyz8mv --temp_dir /emapper_web_jobs/emapper_jobs/user_data/MM_wpfyz8mv --override -m diamond --dmnd_ignore_warnings -i /emapper_web_jobs/emapper_jobs/user_data/MM_wpfyz8mv/queries.fasta --evalue 0.001 --score 60 --pident 40 --query_cover 20 --subject_cover 20 --itype CDS --translate --tax_scope 118882 --target_orthologs all --go_evidence non-electronic --pfam_realign none --report_orthologs --decorate_gff yes --excel > /emapper_web_jobs/emapper_jobs/user_data/MM_wpfyz8mv/emapper.out 2> /emapper_web_jobs/emapper_jobs/user_data/MM_wpfyz8mv/emapper.err
```

最终注释结果中有 7,865 个基因，占总基因 (13,470) 的 58.39%；

注释结果位于 `data/eggnog_mapper_online_result/` 文件夹中；
