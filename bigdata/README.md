##  RefSeq genomes
The reference genomes were downloaded from NCBI with the software [ncbi-genome-download](https://github.com/kblin/ncbi-genome-download) 
and are not included here due to space.
Here is where the downloaded genomes should go, where the `names.txt` files are the output of the ncbi-genome-download runs. 
Notice that the scripts asssume one directory per genus.

If you want to download the Escherichia refseq, you should run something like
```
$ ncbi-genome-download --genus escherichia --assembly-level complete -n --format fasta bacteria | \
 grep GCF | gawk '{print $1}' >> gcf_names.txt
$ ncbi-genome-download --genus escherichia --assembly-level complete --format fasta bacteria
```
Notice that the first command is a dry run (`-n`), which will give you the list of files, that when concatenated to
`gcf_names.txt` will allow for the mapping with the GTDB species names

## GTDB classification
We also dowloaded the classification based on [GTDB](https://gtdb.ecogenomic.org), in particular the 
[metadata from release 89](https://data.ace.uq.edu.au/public/gtdb/data/releases/release89/89.0/bac120_metadata_r89.tsv).
It contains the accession numbers from RefSeq together with the inferred taxonomic classification.
