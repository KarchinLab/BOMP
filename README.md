# BOMP (Burden Or Mutation Position)

BOMP is an association test for case-control study using WES (Whole-Exome Sequencing) data

## Installation

JAVA virtual machine 1.6 or above is required.

## Example

To run BOMP using the example files provided in exmaple/, please put the jar file (bomp.release.jar) and files in example/ in the same directory. Then you should be able to run the following command:

java -Xmx4000M -jar bomp.release.jar -p example.phenotype -v example.var -g example.geno -u example.gene --permutation 1000 --geneSet gene.set

## Parameter setting

Running BOMP without any parameters will list all parameters that can be tuned for BOMP:

java -Xmx4000M -jar bomp.release.jar

Note that the current parameter setting is designed for missense mutations in protein sequences. The algorithm can be applied on DNA sequences of course but the parameters (such as window size) might need to be changed to get better result. For example, the window size of 8 may be too small for DNA sequences.

The defualt number of permutations is 10,000 (which gives the smallest p-value of 1/10001). However, if the test is applied for the whole exome (~20,000 genes), more permutation may be needed to achieve exome-wide significance. For example, 1,000,000 (which gives the smallest p-value of 1/1000001). You can change the number of permutation by --permutation.

## Input file format

### Phenotype file ###

Example: please see example/example.phenotype

Format:
<sample id>,<case-control status (1 for case, 0 for control)>

sample id needs to be unique.

### variant definition file ###

Example: please see eample/example.var

Format:
<variant id>,<gene name>,<gene name (designed for other use in the future)>,<position>,<bioinformatic score>

variant id needs to be unique for each variant. note that column 2 and 3 are the name of the gene where the variant is. the position is for mutation position statistic. since we currently focused on missense mutations, positions are amino acid positions in proteins. A bioinformatics score bewteen 0 (neutral) and 1 (deleterious) indicates its deleteriousness.

### genotype file ###

Example: example/example.geno

Format:
<var 1-sample 1>,<var 2-sample 1>, ... ,<var n-sample 1>
<var 1-sample 2>,<var 2-sample 2>, ... ,<var n-sample 2>
...
<var 1-sample m>,<var 2-sample m>, ... ,<var n-sample m>

Rows represent samples. Columns represent variants. Genotypes are 0 (homozygous ref), 1 (heterozygous), and 2 (homozygous alt). The order of samples needs to be consistent with the order in phenotype file (i.e. k-th row corresponds to the genotypes of the k-th sample in the phenotype file). The order of variants needs to be consistent with the order in variant definition file (i.e. the genotypes in the k-th column corresponds to the variant in the k-th row in the variant definition file).

### gene definition file ###

Example: example/example.gene

Format:
<gene name>,<length>

Make sure that the gene length needs to be greater than the position of any variant occurred in that gene.

### gene set file (optional) ###

Example: example/gene.set

Format:
<set name>\t<description>\t<gene name>[\t<gene name>]

The format follows the text format used in MSigDB (http://www.broadinstitute.org/gsea/msigdb/index.jsp).

