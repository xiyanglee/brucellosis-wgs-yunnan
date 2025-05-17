


```shell
download_eggnog_data.py -P -M -H -d 2 --dbname 'Bacteria' --data_dir /public/data_2/lxy/alto/db/eggnog_mapper_bacteria/ -s
```

```shell
/home/xiyangli/.conda/envs/eggnog-mapper/bin/download_eggnog_data.py:189: SyntaxWarning: invalid escape sequence '\.'
  'outf=$(echo "$file" | sed "s/\.raw_alg\.faa\.gz/\.fa/"); '
Download main annotation database? [y,n] y
Downloading "eggnog.db" at /public/data_2/lxy/alto/db/eggnog_mapper_bacteria...
cd /public/data_2/lxy/alto/db/eggnog_mapper_bacteria && wget -nH --user-agent=Mozilla/5.0 --relative --no-parent --reject "index.html*" --cut-dirs=4 -e robots=off -O eggnog.db.gz http://eggnogdb.embl.de/download/emapperdb-5.0.2/eggnog.db.gz && echo Decompressing... && gunzip eggnog.db.gz
Download taxa database? [y,n] y
Downloading "eggnog.taxa.db" at /public/data_2/lxy/alto/db/eggnog_mapper_bacteria...
cd /public/data_2/lxy/alto/db/eggnog_mapper_bacteria && wget -nH --user-agent=Mozilla/5.0 --relative --no-parent --reject "index.html*" --cut-dirs=4 -e robots=off -O eggnog.taxa.tar.gz http://eggnogdb.embl.de/download/emapperdb-5.0.2/eggnog.taxa.tar.gz && echo Decompressing... && tar -zxf eggnog.taxa.tar.gz && rm eggnog.taxa.tar.gz
Download diamond database (~4GB after decompression)? [y,n] y
Downloading fasta files " at /public/data_2/lxy/alto/db/eggnog_mapper_bacteria...
cd /public/data_2/lxy/alto/db/eggnog_mapper_bacteria && wget -nH --user-agent=Mozilla/5.0 --relative --no-parent --reject "index.html*" --cut-dirs=4 -e robots=off -O eggnog_proteins.dmnd.gz http://eggnogdb.embl.de/download/emapperdb-5.0.2/eggnog_proteins.dmnd.gz && echo Decompressing... && gunzip eggnog_proteins.dmnd.gz
Skipping novel families diamond and annotation databases (or already present). Use -F and -f to force download
Download pfam database (~3GB after decompression)? [y,n] y
Downloading Pfam files " at /public/data_2/lxy/alto/db/eggnog_mapper_bacteria...
cd /public/data_2/lxy/alto/db/eggnog_mapper_bacteria && wget -nH --user-agent=Mozilla/5.0 --relative --no-parent --reject "index.html*" --cut-dirs=4 -e robots=off -O pfam.tar.gz http://eggnogdb.embl.de/download/emapperdb-5.0.2/pfam.tar.gz && echo Decompressing... && tar -zxf pfam.tar.gz && rm pfam.tar.gz
Download MMseqs2 database (~10GB after decompression)? [y,n] y
Downloading MMseqs2 files " at /public/data_2/lxy/alto/db/eggnog_mapper_bacteria...
cd /public/data_2/lxy/alto/db/eggnog_mapper_bacteria && wget -nH --user-agent=Mozilla/5.0 --relative --no-parent --reject "index.html*" --cut-dirs=4 -e robots=off -O mmseqs.tar.gz http://eggnogdb.embl.de/download/emapperdb-5.0.2/mmseqs.tar.gz && echo Decompressing... && tar -zxf mmseqs.tar.gz && rm mmseqs.tar.gz
Download HMMER database of tax ID 2? [y,n] y
Downloading HMMER database of tax ID 2 as "Bacteria" to /public/data_2/lxy/alto/db/eggnog_mapper_bacteria/hmmer/Bacteria
Note that this can take a long time for large taxonomic levels
cd /public/data_2/lxy/alto/db/eggnog_mapper_bacteria/hmmer/Bacteria; echo Downloading HMMs... && wget -N -nH --user-agent=Mozilla/5.0 --relative -r --no-parent --reject "index.html*" --cut-dirs=4 -e robots=off http://eggnog5.embl.de/download/eggnog_5.0/per_tax_level/2//2_hmms.tar.gz && echo Decompressing HMMs... && tar zxf 2_hmms.tar.gz && echo 2/* | xargs mv -t ./ && rm -r 2 && rm 2_hmms.tar.gz; numf=$(find ./ | grep -c ".hmm$"); curr=0; cat /dev/null > Bacteria.hmm_tmp; for file in $(find ./ | grep ".hmm$"); do curr=$((curr+1)); echo "merging HMMs... ${file} (${curr}/${numf})"; cat "${file}" | sed -e "s/.faa.final_tree.fa//" -e "s/.faa.final_tree//" >> Bacteria.hmm_tmp; rm "${file}"; done; mv Bacteria.hmm_tmp Bacteria.hmm; (if [ -f Bacteria.hmm.h3i ]; then rm Bacteria.hmm.h3*; fi) && echo "hmmpress-ing HMMs... " && /home/xiyangli/.conda/envs/eggnog-mapper/bin/hmmpress Bacteria.hmm && echo "generating idmap file... " && cat Bacteria.hmm | grep "^NAME" | sed -e "s/^NAME *//" | awk '{print NR" "$0}' > Bacteria.hmm.idmap && echo "removing single OG hmm files... " && echo ./*hmm | xargs rm;
cd /public/data_2/lxy/alto/db/eggnog_mapper_bacteria/hmmer/Bacteria; echo Downloading FASTAs... && wget -N -nH --user-agent=Mozilla/5.0 --relative -r --no-parent --reject "index.html*" --cut-dirs=4 -e robots=off http://eggnog5.embl.de/download/eggnog_5.0/per_tax_level/2//2_raw_algs.tar && echo Decompressing FASTAs... && tar xf 2_raw_algs.tar && echo 2/* | xargs mv -t ./ && rm -r 2 && rm 2_raw_algs.tar; numf=$(find ./ | grep -c ".faa.gz$"); curr=0; for file in $(find ./ | grep ".faa.gz$"); do curr=$((curr+1)); echo "processing FASTAs...  ${file} (${curr}/${numf})"; outf=$(echo "$file" | sed "s/\.raw_alg\.faa\.gz/\.fa/"); zcat "$file" | awk '/^>/{print; next}{gsub("-", ""); print}' > "$outf" && rm "$file"; done
Finished.
```


```shell
# Download main annotation database
wget -nH --user-agent=Mozilla/5.0 --relative --no-parent --reject "index.html*" --cut-dirs=4 -e robots=off -O eggnog.db.gz http://eggnogdb.embl.de/download/emapperdb-5.0.2/eggnog.db.gz
# Download taxa database
wget -nH --user-agent=Mozilla/5.0 --relative --no-parent --reject "index.html*" --cut-dirs=4 -e robots=off -O eggnog.taxa.tar.gz http://eggnogdb.embl.de/download/emapperdb-5.0.2/eggnog.taxa.tar.gz
# Download diamond database (~4GB after decompression)
wget -nH --user-agent=Mozilla/5.0 --relative --no-parent --reject "index.html*" --cut-dirs=4 -e robots=off -O eggnog_proteins.dmnd.gz http://eggnogdb.embl.de/download/emapperdb-5.0.2/eggnog_proteins.dmnd.gz
# Download pfam database (~3GB after decompression)
wget -nH --user-agent=Mozilla/5.0 --relative --no-parent --reject "index.html*" --cut-dirs=4 -e robots=off -O pfam.tar.gz http://eggnogdb.embl.de/download/emapperdb-5.0.2/pfam.tar.gz
# Download MMseqs2 database (~10GB after decompression)
wget -nH --user-agent=Mozilla/5.0 --relative --no-parent --reject "index.html*" --cut-dirs=4 -e robots=off -O mmseqs.tar.gz http://eggnogdb.embl.de/download/emapperdb-5.0.2/mmseqs.tar.gz
# Download HMMER database of tax ID 2
wget -N -nH --user-agent=Mozilla/5.0 --relative -r --no-parent --reject "index.html*" --cut-dirs=4 -e robots=off http://eggnog5.embl.de/download/eggnog_5.0/per_tax_level/2//2_hmms.tar.gz
wget -N -nH --user-agent=Mozilla/5.0 --relative -r --no-parent --reject "index.html*" --cut-dirs=4 -e robots=off http://eggnog5.embl.de/download/eggnog_5.0/per_tax_level/2//2_raw_algs.tar
```

```shell
cd /public/data_2/lxy/alto/db/eggnog_mapper_bacteria

# Download main annotation database
gunzip eggnog.db.gz

# Download taxa database
tar -zxf eggnog.taxa.tar.gz && rm eggnog.taxa.tar.gz

# Download diamond database (~4GB after decompression)
gunzip eggnog_proteins.dmnd.gz

# Download pfam database (~3GB after decompression)
tar -zxf pfam.tar.gz && rm pfam.tar.gz

# Download MMseqs2 database (~10GB after decompression)
tar -zxf mmseqs.tar.gz && rm mmseqs.tar.gz

# Download HMMER database of tax ID 2
mkdir -p ./hmmer/Bacteria
mv 2_hmms.tar.gz 2_raw_algs.tar ./hmmer/Bacteria/
cd ./hmmer/Bacteria

tar zxf 2_hmms.tar.gz && echo 2/* | xargs mv -t ./ && rm -r 2 && rm 2_hmms.tar.gz; numf=$(find ./ | grep -c ".hmm$"); curr=0; cat /dev/null > Bacteria.hmm_tmp; for file in $(find ./ | grep ".hmm$"); do curr=$((curr+1)); echo "merging HMMs... ${file} (${curr}/${numf})"; cat "${file}" | sed -e "s/.faa.final_tree.fa//" -e "s/.faa.final_tree//" >> Bacteria.hmm_tmp; rm "${file}"; done; mv Bacteria.hmm_tmp Bacteria.hmm; (if [ -f Bacteria.hmm.h3i ]; then rm Bacteria.hmm.h3*; fi) && echo "hmmpress-ing HMMs... " && /home/xiyangli/.conda/envs/eggnog-mapper/bin/hmmpress Bacteria.hmm && echo "generating idmap file... " && cat Bacteria.hmm | grep "^NAME" | sed -e "s/^NAME *//" | awk '{print NR" "$0}' > Bacteria.hmm.idmap && echo "removing single OG hmm files... " && echo ./*hmm | xargs rm;

tar xf 2_raw_algs.tar && echo 2/* | xargs mv -t ./ && rm -r 2 && rm 2_raw_algs.tar; numf=$(find ./ | grep -c ".faa.gz$"); curr=0; for file in $(find ./ | grep ".faa.gz$"); do curr=$((curr+1)); echo "processing FASTAs...  ${file} (${curr}/${numf})"; outf=$(echo "$file" | sed "s/\.raw_alg\.faa\.gz/\.fa/"); zcat "$file" | awk '/^>/{print; next}{gsub("-", ""); print}' > "$outf" && rm "$file"; done
```


