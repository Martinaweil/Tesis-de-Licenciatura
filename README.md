## Caracterización funcional de los 'genes seminales' en especies cactófilas del género *Drosophila*
Este repositorio tiene todos los archivos y scripts utilizados en el proceso de análisis de datos de mi tesis de licenciatura.

**1-Instalaciones**

```bash
curl -O https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
sh Miniconda3-latest-Linux-x86_64.sh   #instalacion de bioconda para linux

conda config --add channels bioconda     #set up channels
conda config --add channels conda-forge   #necesario para descargar paquetes

#Instalación de paquetes usados

conda install -c bioconda trinity #trinity
conda install -c bioconda trinotate  #trinotate
conda install -c bioconda bioconductor-goseq  #goseq
conda install -c bioconda orthofinder  #orthofinder
```

**2-Trinotate**

Trinotate se utiliza para anotar funcionalmente el transcriptoma deseado, en este caso el transcriptoma de D *buzzatii*. 

```bash
exec bash 
Build_Trinotate_Boilerplate_SQLite_db.pl  Trinotate #descarga bases de datos y arma base de datos SQL de trinotate
makeblastdb -in uniprot_sprot.pep -dbtype prot #crea base de datos en formato para BLAST a partir de base de datos uniprot
gunzip Pfam-A.hmm.gz   #preparo Pfam para el uso con HMMER
hmmpress Pfam-A.hmm

#BLAST
#blastx: toma transcriptos como input, los traduce y busca en bd proteica
blastx -query transcripts.fasta -db uniprot_sprot.pep -num_threads 4 -max_target_seqs 1 -outfmt 6 -evalue 1e-5 > blastx.outfmt6 

#blastp: toma proteínas como input y busca en bd proteica
blastp -query Dbor2/transcripts.fasta.transdecoder.pep -db uniprot_sprot.pep -num_threads 4 -max_target_seqs 1 -outfmt 6 -evalue 1e-5 > Dbor2/blastp.outfmt6 

#HMMER
#toma proteínas como input y busca en base de datos de dominios proteicos. 
hmmscan --cpu 4 --domtblout Dbor2/TrinotatePFAM.out Pfam-A.hmm  Dbor2/transcripts.fasta.transdecoder.pep > pfam.log 

#cargo la base de datos “Trinotate.sqlite” con el mapa de genes y transcriptos que me da trinity, mis transcriptos, y mi output de transdecoder
Trinotate Dkoe2/Trinotate_V3.sqlite init --gene_trans_map Dkoe2/transcripts.gene_trans_map --transcript_fasta Dkoe2/transcripts.fasta --transdecoder_pep Dkoe2/transcripts.fasta.transdecoder.pep 

#cargo los outputs de blastp y blastx a la base de datos
Trinotate Dkoe2/Trinotate_V3.sqlite LOAD_swissprot_blastp Dkoe2/blastp.outfmt6
Trinotate Dkoe2/Trinotate_V3.sqlite LOAD_swissprot_blastx Dkoe2/blastx.outfmt6

#cargo output de hmmer
Trinotate Dkoe2/Trinotate_V3.sqlite LOAD_pfam Dkoe2/TrinotatePFAM.out

#genero output de trinotate
Trinotate Dkoe2/Trinotate_V3.sqlite report [opts] > Dkoe2/trinotate_annotation_report.xls

```

**3-GOseq**

Lo que quiero hacer es:

```bash
anaconda3/opt/trinity-2.1.1/Analysis/DifferentialExpression/run_GOseq.pl \
                       --factor_labeling  Dkoe2/factor_labeling.txt \
                       --GO_assignments Dkoe2/go_annotations.txt \
                       --lengths Dkoe2/Trinity.gene_lengths.txt \
                       --background Dkoe2/genes_ids.txt
```

Por lo tanto primero necesito generar esos documentos.

1. go_annotations.txt

```bash
Exec bash

#Extraigo los GO terms de mi output de trinotate

extract_GO_assignments_from_Trinotate_xls.pl\
                         --Trinotate_xls Dkoe2/trinotate_annotation_report.xls \
                         -G --include_ancestral_terms \
                         > Dkoe2/go_annotations.txt
```

2. Trinity.gene_lengths.

```bash
exec bash
fasta_seq_length.pl  Dkoe2/transcripts.fasta > Dkoe2/Trinity.fasta.seq_lens
#y a partir de este archivo:
# esto lo tengo que hacer en python 2.7!
conda activate python2
TPM_weighted_gene_length.py  \
         --gene_trans_map Dkoe2/transcripts.gene_trans_map \
         --trans_lengths Dkoe2/Trinity.fasta.seq_lens \
         --TPM_matrix Dkoe2/RSEM.trans.TMM.EXPR.matrix > Dkoe2/Trinity.gene_lengths.txt
```

3. Factor_labeling.txt:

Es un archivo en formato factor (tab) gene_id, donde “factor” describe un grupo de genes. En nuestro caso control y glándulas accesorias



**4-Orthofinder**

Orthofinder se utilizó para ampliar la lista de ACPs de D *buzzatii* y para analizar si los ACPs de D *virilis* y D *buzzatii* son ortólogos. 

```R
#En R preparo las listas de acps que necesito

#leo las listas de genes de ACPs que estan en factor labeling
Dant_acps <- read.delim("~/tesis-gene-annotation/Especies/Dant_basedatos_drosophila/factor_labeling.txt", header=FALSE)
Dbuz_acps <- read.delim("~/tesis-gene-annotation/Especies/Dbuz_basedatos_drosophila/factor_labeling.txt", header=FALSE)
Dbor_acps <- read.delim("~/tesis-gene-annotation/Especies/Dbor_basedatos_drosophila/factor_labeling.txt", header=FALSE)
Dkoe_acps <- read.delim("~/tesis-gene-annotation/Especies/Dkoe_basedatos_drosophila/factor_labeling.txt", header=FALSE)

#leo el documento donde se mapea la/s proteina/s correspondiente/s
Dant <- read.delim("~/tesis-gene-annotation/Especies/Dant_basedatos_drosophila/trinotate_annotation_report.xls", header=TRUE)
Dbuz <- read.delim("~/tesis-gene-annotation/Especies/Dbuz_basedatos_drosophila/trinotate_annotation_report.xls", header=TRUE)
Dbor <- read.delim("~/tesis-gene-annotation/Especies/Dbor_basedatos_drosophila/trinotate_annotation_report.xls", header=TRUE)
Dkoe <- read.delim("~/tesis-gene-annotation/Especies/Dkoe_basedatos_drosophila/trinotate_annotation_report.xls", header=TRUE)

#genero tabla de IDs de proteinas provenientes de genes que codifican ACPs
pequeña_tabla_Dant<-subset(Dant, gene_id %in% Dant_acps$V2)
pequeña_tabla_Dbuz<-subset(Dbuz, gene_id %in% Dbuz_acps$V2)
pequeña_tabla_Dbor<-subset(Dbor, gene_id %in% Dbor_acps$V2)
pequeña_tabla_Dkoe<-subset(Dkoe, gene_id %in% Dkoe_acps$V2)

DBUZ<-data.frame(pequeña_tabla_Dbuz[,2])
DBUZ<-paste0("Dbuz_",DBUZ$pequeña_tabla_Dbuz...2.)
DBUZ<-data.frame(DBUZ)
setwd("~/tesis-gene-annotation/orthofinder/selected_proteins/selected_proteins/OrthoFinder/Results_Aug04_1")
write.csv(DBUZ,file = "Dbuz_acps.txt",quote = FALSE,row.names = FALSE,col.names = FALSE)

```

En una carpeta pongo los archivos FASTA de proteinas para cada especie y corro orthofinder:

```bash
orthofinder -f *nombredelacarpeta*/
```

Esto lo realicé dos veces:

- Por un lado para el cluster Buzzatii y D *mojavensis* juntos. A partir de esta corrida determine el nuevo grupo de ACPs de D *buzzatii* bajo la condición que sea una proteína perteneciente a un ortogrupo donde hay al menos dos ACPs hipotéticos. 

Para filtrar ortogrupos uso el script og_filtering.R usando el siguiente comando:

```bash
Rscript --vanilla  og_filtering.R -g Orthogroups/Orthogroups.txt -l *archivo_lista_ACPs*.txt  -n *filter_name*
```

Cuando quiero filtrar por mas grupos tengo que cambiarle el nombre a la carpeta "filtered_orthogroups" para que no borre todo y la vuelva a crear. Asimismo le debo agregar al comando anterior "-f filtered_groups.rds"

Luego quiero procesar el output de og_filtering.r

```R
#Procesa el output del filtered_orthogroups buscando todos los orthogroups que tengan almenos dos ACPS de especies distintas

data<-read.table("filtered_orthogroups/annotated_groups.txt",header =TRUE)
colnames(data) <- c("Orthogroup","Dant", "Dbuz", "Dbor", "Dkoe")
apply(data[,2:4], 1, any)
data2 <- which(rowSums(data[,2:5])>1)
data2<-data[data2, ]
#write.table(data[data2, ],file = "annotated_groups_FILTRADO_MIN_DOS_TRUE.txt",quote = FALSE,row.names = FALSE)


#ahora quiero filtrar el orthogroups y tomar solo aquellos que aparezcan en mi lista de filtrados con al menso dos true

ortho<-read.delim("Orthogroups/orthogroups_Dbuz")
filtrado<-ortho[ortho$Orthogroup %in% data2$Orthogroup,]

write.table(filtrado,file = "nuevosacps_Dbuz.txt",quote = FALSE,row.names = FALSE)


```

A partir de esta corrida de Orthofinder genero también una nueva lista de genes de AG que voy a necesitar para correr GOseq:

```R
#luego de orthofinder
nuevosacps_Dbuz <- read.csv("~/Dbuz_nuevosAG/nuevosacps_DbuzESTRICTO.txt", sep="",header = FALSE)
#creo mi lista de genes AG para usar en goseq
nuevos_AG<-subset(Dbuz, transcript_id %in% nuevosacps_Dbuz$V1)
nuevos_AG<-(unique(nuevos_AG$gene_id))
V2<-c(rep("AG",47))
fact<-cbind(V2,nuevos_AG)
write.table(fact,file = "AG_estricto.txt",quote=FALSE, row.names = FALSE, col.names = FALSE,sep = "\t" )
 
```

- Por otro lado para D buzzatii con D *mojavensis*, D *melanogaster* y D *virilis* para ver si los ACPs de D *virilis* y D *buzzatii* son ortólogos. 

Para crear diagramas de Venn:

```R
annotated_groups <- read.delim("~/tesis-gene-annotation/orthofinder/Orthofinder_2_DbuzDmelDmojDvir/OrthoFinder/Results_Sep03/filtered_orthogroups/annotated_groups.txt")
#queremos hacer dos vectores, uno para cada especie, donde figuren los ortogrupos que tengan TRUE para luego realizar el diagrama de Venn.

Dvir<-annotated_groups[ which(annotated_groups$Dvir==TRUE),"ogID"]
Dbuz<-annotated_groups[ which(annotated_groups$Dbuz==TRUE),"ogID"]

venn.diagram(
  x=list(Dvir,Dbuz),
  category.names = c("D. virilis","D buzzatii"),
  filename = "ortogrupos_vir+buz.png"
)

```

