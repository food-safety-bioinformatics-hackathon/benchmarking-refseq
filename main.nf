#!/usr/bin/env nextflow
// nextflow.preview.dsl = 2
params.ref = "/qib/platforms/Informatics/hackathon2019/data/NZ_CP009072.1.fasta"
params.length_reads = 150
params.seed = 1000
params.output = "results"
params.rank = false
params.score = false

refs =  Channel.fromPath(params.ref)

process simulation {
publishDir "${params.output}/simulation", mode: "copy"
tag {seed}
input:
file(ref) from refs
each sim from 1..10
output:
set seed, file("sim_${name}_${seed}_1.fq"), file("sim_${name}_${seed}_2.fq") into (sim_ch_shovill, sim_ch_refrank, sim_ch_refseq)

script:
seed = params.seed + sim
name = ref.getBaseName()
"""
wgsim -1 ${params.length_reads} -2 ${params.length_reads} -S ${seed} ${ref} sim_${name}_${seed}_1.fq sim_${name}_${seed}_2.fq
"""
}

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
refseq_masher matches --top-n-results 50 --output "refseq_masher_assembly_${seed}.tsv" output/contigs.fa
"""
}

process refseq_masher_reads {
publishDir "${params.output}/refseq_masher_reads", mode: "move"

cpus 8

input:
set seed, file(r1), file(r2) from sim_ch_refseq

output:
file("refseq_masher_reads_${seed}.tsv") into masher_reads_ch

script:

"""
refseq_masher contains --top-n-results 50 --parallelism ${task.cpus} --output "refseq_masher_reads_${seed}.tsv" $r1 $r2
"""
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

if (params.score) {
process score {
input:
file("refseq_masher_reads_${seed}.tsv") from 
file("refseq_masher_assembly_${seed}.tsv")
output:

script:

"""
#!/usr/bin/python
with open("${paraarams.ref}.out", "w") as out_fh:
  with open("refseq_masher_reads_${seed}.tsv", "r") as in_fh:
    lines = fh.readlines()
"""
  }
}

