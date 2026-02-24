import os
import glob
from snakemake.io import glob_wildcards
from snakemake.io import expand

DATA_DIR=config["Data_dir"]
REF_PATH=config["Ref_path"]
META_PATH=config["Meta_path"]

PATTERN = "{batch}/{prefix}_MiSeq-{suffix}/{sample}.fastq.gz"
samples_path = os.path.join(DATA_DIR, PATTERN)
batches, prefixes, suffixes, samples = glob_wildcards(samples_path)

run_map = {
  (b, s): f"{p}_MiSeq-{sx}" 
  for b, p, sx, s in zip(batches, prefixes, suffixes, samples)}
PAIRS = sorted(run_map.keys())

#targets = [f"sports_final/{b}/1_{s}/{s}_result/{s}_output.txt" for b, s in PAIRS]

rule all:
    input:
        "out_data/01_data_processing.html",
        "out_data/02_overview_plots.html",
        "out_data/03_de_analysis.html"

 
def input_path(wc):
    run = run_map[(wc.batch, wc.sample)]
    return os.path.join(DATA_DIR, wc.batch, run, f"{wc.sample}.fastq.gz")        

# For the Adapter-seq see: https://support-docs.illumina.com/SHARE/AdapterSequences/Content/SHARE/AdapterSeq/TruSeq/TruSeqSmallRNA.htm
rule trimming:
    input:
      input_path
    output:
      trimmed_fastq="trimmed_final/{batch}/{sample}.fastq"
    params:
      adapter_seq="TGGAATTCTCGGGTGCCAAGG",
      out_dir="trimmed_final"
    threads: 1
    shell:
      """
        mkdir -p {params.out_dir}
        cutadapt -a "{params.adapter_seq}" -o "{params.out_dir}/{wildcards.batch}/{wildcards.sample}.fastq" {input} --maximum-length 50 --minimum-length 15
      """

rule reads_cleaning:
    input:
      trimmed_fastq=rules.trimming.output.trimmed_fastq
    output:
      filtered_fastq="filtered_final/{batch}/{sample}.fastq"
    params:
      trimmed_dir=rules.trimming.params.out_dir,
      out_dir="filtered_final"
    threads: 1
    script: "src/readCleaning.R"
 
rule sports_align:
  input:
    filtered_fastq=rules.reads_cleaning.output.filtered_fastq
  output:
    aligned_sample="sports_final/{batch}/1_{sample}/{sample}_result/{sample}_output.txt"
  params:
    ref_path=REF_PATH,
    out_dir="sports_final"
  threads: 4
  log: "logs/sports_{batch}_{sample}.txt"
  shell: 
    """
      exec 2>>{log} 
      mkdir -p {params.out_dir}
      sports.pl -i {input.filtered_fastq} -p 8 -o {params.out_dir}/{wildcards.batch} -M 2 \
      -g {params.ref_path}/UCSC/mm10/Sequence/BowtieIndex/genome \
      -m {params.ref_path}/miRBase_21/miRBase_21-mmu \
      -r {params.ref_path}/rRNAdb/mouse_rRNA \
      -t {params.ref_path}/GtRNAdb/mm10-tRNAs \
      -e {params.ref_path}/Ensembl/Mus_musculus.GRCm38.ncrna \
      -f {params.ref_path}/Rfam_12.3/Rfam-12.3-mouse \
      -w {params.ref_path}/piRBase/piR_mouse \
      -L 50
    """

def sports_all_samples(wc):
    return [
        f"sports_final/{b}/1_{s}/{s}_result/{s}_output.txt"
        for (b, s) in PAIRS
        ]

rule construct_se:
  input:
    sports_all_samples
  output:
    output_html="out_data/01_data_processing.html",
    se="out_data/01_sports_se.rds"
  params:
    out_dir="out_data",
    meta_path=META_PATH,
    sports_dir=rules.sports_align.params.out_dir
  threads: 1
  log: "logs/construct_se.log"
  shell:
    """
      exec 2>>{log} 
      mkdir -p {params.out_dir}
      Rscript -e 'rmarkdown::render("src/01_data_processing.Rmd", 
                                    "html_document", 
                                     output_file="../{output.output_html}",
                                     params=list(sports_dir="../{params.sports_dir}",
                                                 out_dir="../{params.out_dir}",
                                                 meta_path="../{params.meta_path}"))'
    """


#TODO: check filtering reads stats first!!    
rule quality_plots:
  input:
    se=rules.construct_se.output.se
  output:
    se="out_data/02_sports_se.rds",
    output_html="out_data/02_overview_plots.html"
  params:
    out_data_dir=rules.construct_se.params.out_dir,
    out_plot_dir="plots"
  threads: 1
  log: "logs/quality_plots.log"
  shell:
    """
      mkdir -p {params.out_plot_dir}
      exec 2>>{log} 
      Rscript -e 'rmarkdown::render("src/02_overview_plots.Rmd", 
                                    "html_document", 
                                     output_file="../{output.output_html}",
                                     params=list(se_dir="../{input.se}",
                                                 out_data_dir="../{params.out_data_dir}",
                                                 out_plots_dir="../{params.out_plot_dir}"))'
    
    """

rule differential_testing: 
  input:
    se=rules.quality_plots.output.se
  output:
    output_html="out_data/03_de_analysis.html"
  params:
    out_data_dir=rules.construct_se.params.out_dir,
    out_plot_dir=rules.quality_plots.params.out_plot_dir,
  threads: 1
  log: "logs/differential_testing.log"
  shell:
    """
      exec 2>>{log} 
      Rscript -e 'rmarkdown::render("src/03_de_analysis.Rmd", 
                                    "html_document", 
                                     output_file="../{output.output_html}",
                                     params=list(se_dir="../{input.se}",
                                                 out_data_dir="../{params.out_data_dir}",
                                                 out_plots_dir="../{params.out_plot_dir}"))'
    """
