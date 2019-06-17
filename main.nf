#!/usr/bin/env nextflow
// nextflow.preview.dsl = 2
params.refdir = "/media/deolivl/QIB_deolivl/bigdata/noplasmids" 
params.length_reads = 150
params.seed = 1003
params.output = "results"
params.rank = false
params.assembly = false
params.sim = 10
params.besthits = 50
params.wgsim_options = "-N 100000"

reflist =  Channel.fromPath(params.refdir + "/*.gz")
ref = reflist.randomSample(1)

process simulation {
  publishDir "${params.output}/simulation", mode: "copy"
  tag {seed}
  input:
    file ref 
    each sim from 1..params.sim
  output:
    set seed, file("sim_${name}_${seed}_1.fq"), file("sim_${name}_${seed}_2.fq") into (sim_ch_shovill, sim_ch_refrank, sim_ch_refseq, sim_ch_get_cov)
    file ref into sim_ref_ch
  script:
    seed = params.seed + sim
    name = ref.getBaseName()
    """
    if [ ! -f ref.fas ]; then  ## must be unzipped
      zcat ${ref} > ref.fas
    fi
    wgsim $params.wgsim_options -1 ${params.length_reads} -2 ${params.length_reads} -S ${seed} ref.fas sim_${name}_${seed}_1.fq sim_${name}_${seed}_2.fq
    """
}

if (params.assembly) {
  process shovill {
  publishDir "${params.output}/shovill", saveAS: {"sim_${seed}"}, mode: "copy"
  cpus 8
  tag {seed}
  input:
    set seed, file(r1), file(r2) from sim_ch_shovill
  output:
    set seed, file("output") into shovill_ch
  script:
    """
    shovill --outdir output --R1 $r1 --R2 $r1
    """
  }

  process refseq_masher_assembly {
  publishDir "${params.output}/refseq_masher_assembly"
  tag {seed}
  input:
    set seed, file(output) from shovill_ch
  output:
    file("refseq_masher_assembly_${seed}.tsv") into masher_assembly_ch
  script:
    """
    refseq_masher matches --top-n-results 200 --output "refseq_masher_assembly_${seed}.tsv" output/contigs.fa
    """
  }
}

process refseq_masher_reads {
  publishDir "${params.output}/refseq_masher_reads", mode: "copy"
  cpus 4
  input:
    set seed, file(r1), file(r2) from sim_ch_refseq
  output:
    file("refseq_masher_reads_${name}.tsv") into masher_reads_ch
  script:
    name = r1.getBaseName()
    """
    refseq_masher contains --top-n-results 500 --parallelism ${task.cpus} --output "refseq_masher_reads_${name}.tsv" $r1 $r2
    """
}

process find_GCF {
	publishDir "${params.output}/tophits", mode: "copy" 
  input:
    file tsv from masher_reads_ch
  output:
    file("${name}_hits.tsv") into find_gcf_ch
  
  script:
    name = tsv.getBaseName()

    """
    #!/usr/bin/env python
    import re
    with open("${tsv}","r") as ifh:
      lines = ifh.readlines()
      headers = lines[0].split("\\t")
      assembly_accession_index = headers.index("assembly_accession")
      hits = []
      for line in lines[1:]:
          asm = line.split("\\t")[assembly_accession_index]
          if (re.match("^GCF", asm)):
              hits.append(line.split("\\t")[assembly_accession_index])
      with open("${name}_hits.tsv", "w") as ofh:
          for hit in hits:
              ofh.write("{}\\n".format(hit))
    """
}

process get_coverage {
  publishDir "${params.output}/tophits", mode: "copy"
  tag {name}

  input:
    set seed, file(r1), file(r2) from sim_ch_get_cov
    params.besthits
    file(tophits) from find_gcf_ch
    file(ref) from sim_ref_ch
    val params.refdir
  output:
    file("final_scores_${name}_${seed}.csv")
    file("debug_${seed}.csv")
  shell:
    name = ref.getBaseName()
    '''
    typeset -i count=0
    for gcf in `cat !{tophits}`; do
      if [ $count -gt !{params.besthits} ]; then break;
      else
        HITREF=`\\ls !{params.refdir}/${gcf}*`
        if test x"${HITREF}" != x""; then
          echo ${HITREF} >> debug_!{seed}.csv
          count=$count+1
          bwa index ${HITREF}
          bwa mem ${HITREF} !{r1} !{r2} | samtools view -S -b - | samtools sort > sort.bam
          bamcov -H sort.bam >> scores.csv
        fi  # if HITREF exists
      fi    # if count == 10
    done
    bwa index !{ref}
    bwa mem !{ref} !{r1} !{r2} | samtools view -S -b - | samtools sort > sort.bam
    bamcov -H sort.bam >> scores.csv
    score.py scores.csv > final_scores_!{name}_!{seed}.csv
    '''
}

if (params.rank) {
  process ref_rank {
  input:
    set seed, file(r1), file(r2) from sim_ch_refrank
  output:

  script:

    """
    refrank.py --ref *fna --fastq sim_${sim}_1.fq sim_${sim}_2.fq
    """
	}
}

