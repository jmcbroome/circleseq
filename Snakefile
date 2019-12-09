rule all:
    input:
        "{sample}_{reference}_errors.png"
rule bwa_map:
    input:
        "input/{sample}_R1.fa",
        "input/{sample}_R2.fa",
        "references/{reference}.fa"
    output:
        "{sample}_{reference}_initial_mapping.sam"
    benchmark:
        "{sample}_{reference}_imap_benchmark.tsv"
    threads: 20
    shell: 
        "bwa mem -t {threads} {input[2]} {input[0]} {input[1]} > {output}"
rule process_circle_alignments:
    input:
        "{sample}_{reference}_initial_mapping.sam"
    output:
        "{sample}_{reference}_split.fa",
        "{sample}_{reference}_regions.bed"
    benchmark:
        "{sample}_{reference}_process_benchmark.tsv"
    shell:
        "python3 process_circle_alignments.py -f {output[0]} -b {output[1]} -d 10 {input}"
rule makefasta:
    input:
        "{sample}_{reference}_regions.bed",
        "references/{reference}.fa"
    output:
        "{sample}_{reference}_reference.fa"
    shell:
        "bedtools getfasta -fi {input[1]} -bed {input[0]} -fo {output} -name"
rule make_brpileup:
    input:
        "{sample}_{reference}_split.fa",
        "{sample}_{reference}_regions.bed",
        "{sample}_{reference}_reference.fa"
    output:
        "{sample}_{reference}_consensus.sam"
    benchmark:
        "{sample}_{reference}_refpileup_benchmark.tsv"
    threads: 20
    shell:
        "python3 make_brpileup.py -t {threads} -c {output} {input[1]} {input[2]} {input[0]}"
rule index:
    input:
        "references/{reference}.fa"
    output:
        "references/{reference}.fa.fai"
    shell:
        "samtools faidx references/{wildcards.reference}.fa"
rule samprocess:
    input:
        "{sample}_{reference}_consensus.sam",
        "references/{reference}.fa.fai"
    output:
        "{sample}_{reference}_sorted.bam"
    shell:
        "samtools view -ht {input[1]} {input[0]} | samtools view -Sb | samtools sort > {output}"
rule mpileup:
    input:
        "{sample}_{reference}_sorted.bam"
    output:
        "{sample}_{reference}_variants.txt"
    shell:
        "samtools mpileup -f references/{wildcards.reference}.fa {input} > {output}"
rule count_coverage:
    input:
        "{sample}_{reference}_sorted.bam",
        "references/{reference}.fa"
    output:
        "{sample}_{reference}_coveredbases.fa"
    run:
        shell("bedtools genomecov -bg -ibam {input[0]} > {wildcards.sample}_{wildcards.reference}_coverage.bedgraph")
        shell("bedtools getfasta -fi {input[1]} -bed {wildcards.sample}_{wildcards.reference}_coverage.bedgraph -fo {output}") 
rule count_errors:
    input:
        "{sample}_{reference}_variants.txt",
        "{sample}_{reference}_coveredbases.fa"
    output:
        "{sample}_{reference}_errors.txt"
    shell:
        "python3 count_errors.py -t 2 -e {input[0]} -s {input[1]} > {output}"
rule graph_errors:
    input:
        "{sample}_{reference}_errors.txt"
    output:
        "{sample}_{reference}_errors.png"
    shell:
        "python3 graph_errors.py {input}"