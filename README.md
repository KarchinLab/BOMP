# BOMP (Burden Or Mutation Position)

BOMP is an association test for case-control study using WES (Whole-Exome Sequencing) data

## Installation

JAVA virtual machine 1.6 or above is required.

## Example

Put **bomp.release.jar** and **files in example/** in the same directory.

Run BOMP using the following command:

java -Xmx4000M -jar bomp.release.jar -p example.phenotype -v example.var -g example.geno -u example.gene --permutation 1000 --geneSet gene.set

The result will be saved in bomp.out.

## Parameter Setting

Running BOMP without any parameters will list all parameters that can be tuned for BOMP:

**java -Xmx4000M -jar bomp.release.jar**

Note that the current parameter setting is designed for missense mutations in protein sequences. The algorithm can be applied on DNA sequences of course but the parameters (such as window size) might need to be changed to get better result. For example, the window size of 8 may be too small for DNA sequences.

The defualt number of permutations is 10,000 (which gives the smallest p-value of 1/10001). However, if the test is applied for the whole exome (~20,000 genes), more permutation may be needed to achieve exome-wide significance. For example, 1,000,000 (which gives the smallest p-value of 1/1000001). You can change the number of permutation by --permutation.

## File Format

### Phenotype file ###

Example: please see example/example.phenotype

Format:

**sample_id,case-control_status (1 for case, 0 for control)**

sample id needs to be unique.

### Variant definition file ###

Example: please see eample/example.var

Format:

**variant_id,gene_name,gene_name (designed for other use in the future),position,bioinformatic_score**

variant id needs to be unique for each variant. note that column 2 and 3 are the name of the gene where the variant is. the position is for mutation position statistic. since we currently focused on missense mutations, positions are amino acid positions in proteins. A bioinformatics score bewteen 0 (neutral) and 1 (deleterious) indicates its deleteriousness.

### Genotype file ###

Example: example/example.geno

Format:

**var1_for_sample1,var2_for_sample1, ... ,varN_for_sample1**

**var1_for_sample2,var2_for_sample2, ... ,varN_for_sample2**

**...**

**var1_for_sampleM,var2_for_sampleM, ... ,varN_for_sampleM**

Rows represent samples. Columns represent variants. Genotypes are 0 (homozygous ref), 1 (heterozygous), and 2 (homozygous alt). The order of samples needs to be consistent with the order in phenotype file (i.e. k-th row corresponds to the genotypes of the k-th sample in the phenotype file). The order of variants needs to be consistent with the order in variant definition file (i.e. the genotypes in the k-th column corresponds to the variant in the k-th row in the variant definition file).

### Gene definition file ###

Example: example/example.gene

Format:

**gene_name,length**

Make sure that the gene length needs to be greater than the position of any variant occurred in that gene.

### Gene set file (optional) ###

Example: example/gene.set

Format:

**gene_set_name \t description \t gene_name \[\t gene name\]

The format follows the text format used in MSigDB (http://www.broadinstitute.org/gsea/msigdb/index.jsp).

