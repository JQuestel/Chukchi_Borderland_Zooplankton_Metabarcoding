set.current(processors=10)
make.file(type=fastq, prefix=V4)
make.contigs(file=current, insert=30) 
summary.seqs(fasta=current)
screen.seqs(fasta=current, group=current, summary=current, minlength=300, maxambig=0)
summary.seqs(fasta=current)
unique.seqs(fasta=current)
count.seqs(name=current, group=current)
summary.seqs(fasta=current, name=current)
align.seqs(fasta=current, reference=132_for_mothur.V4pcr.fasta, flip=t)
summary.seqs(fasta=current, count=current)
get.current()
screen.seqs(fasta=current, count=current, summary=current, start=24, end=8723, maxhomop=8)
filter.seqs(fasta=current, vertical=T)
unique.seqs(fasta=current, count=current)
pre.cluster(fasta=current, count=current, diffs=2, method=unoise)
chimera.vsearch(fasta=current, count=current, dereplicate=t)
remove.seqs(fasta=current, accnos=current, dups=f)
summary.seqs(fasta=current, count=current)
get.current()
classify.seqs(fasta=current, count=current, reference=132_for_mothur.V4pcr.ng.fasta, taxonomy=taxmap_ssu_132_Ready.V4pcr.txt, cutoff=80)
make.shared(count=current, label=asv)
classify.otu(list=current, count=current, taxonomy=current, label=asv)
get.current()
#remove non-mesozooplankton taxa
remove.lineage(fasta=current, count=current, taxonomy=current, taxon=Eukaryota;Archaeplastida;-Eukaryota;Centrohelida;-Eukaryota;Eukaryota_unclassified;-Eukaryota;Haptophyta;-Eukaryota;Incertae_Sedis;-Eukaryota;Opisthokonta;Holozoa;Choanoflagellida;-Eukaryota;Opisthokonta;Holozoa;Filasterea;-Eukaryota;Opisthokonta;Holozoa;Holozoa_unclassified;-Eukaryota;Opisthokonta;Holozoa;Metazoa_|Animalia|;Eumetazoa;Bilateria;Bilateria_unclassified;-Eukaryota;Opisthokonta;Holozoa;Metazoa_|Animalia|;Eumetazoa;Bilateria;Vertebrata;Tetrapoda;-Eukaryota;Opisthokonta;Holozoa;Metazoa_|Animalia|;Eumetazoa;Eumetazoa_unclassified;-Eukaryota;Opisthokonta;Holozoa;Metazoa_|Animalia|;Metazoa_|Animalia|_unclassified;-Eukaryota;Opisthokonta;Opisthokonta_unclassified;-Eukaryota;SAR;)
summary.tax(taxonomy=current, count=current)
summary.seqs(fasta=current, count=current)
#need version 1.44 or newer 
make.shared(count=current, label=asv)
#re-classify ASVs to get new cons.taxonomy files
classify.otu(list=current, count=current, taxonomy=current, label=asv)
count.groups(shared=current)
#removes rare OTUs; global singletons occurring <2 
remove.rare(list=current, count=current, nseqs=2, label=asv)
make.shared(count=current, label=asv)
get.oturep(fasta=current, count=current, list=current, method=abundance)
summary.seqs(fasta=current, count=current)
classify.otu(taxonomy=current, count=current, list=current)

#alpha diversities
summary.single(shared=current, calc=nseqs-sobs-coverage-shannon-shannoneven-invsimpson, label=asv)
#beta diversity - Bray Curtis
dist.shared(shared=current, calc=braycurtis, subsample=T, label=asv)

#####Copepod IDs using ArCop
classify.seqs(fasta=current, count=current, reference=18SArcticCopepodaWorkingAligned.V4pcr.ng.fasta, taxonomy=TaxMap_Mothur_18SArcticCopepoda.V4pcr.txt, cutoff=80)
make.shared(count=current, label=asv)
classify.otu(taxonomy=current, count=current, list=current)

