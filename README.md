# Conpair
Conpair: concordance and contamination estimator for tumor–normal pairs

Conpair is a fast and robust method dedicated for human tumor-normal studies to perform concordance verification (= samples coming from the same individual), as well as cross-individual contamination level estimation in whole-genome and whole-exome sequencing experiments. Importantly, our method of estimating contamination in the tumor samples is not affected by copy number changes and is able to detect contamination levels as low as 0.1%.

* Version: 0.2
* Author: Ewa A Bergmann
* Contact: ewa.a.bergmann@gmail.com

**Required input files:** two bam files (tumor, normal)

**Required software:** GATK 2.3 or later, python 2.7 or higher, scipy, numpy, java

**Required data:** Human genome file (GRCh37 or GRCh38) 

GRCh37:

The fasta file can be downloaded from: ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37.fasta.gz  
In order to be able use the fasta file as a reference 2 additional files are required:
`human_g1k_v37.dict`, `human_g1k_v37.fa.fai`  
To create these files please follow: http://gatkforums.broadinstitute.org/gatk/discussion/1601/how-can-i-prepare-a-fasta-file-to-use-as-reference

GRCh38:

The fasta file can be downloaded from: http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
In order to be able use the fasta file as a reference 2 additional files are required:
`hg38.dict`, `hg38.fa.fai`
To create these files please follow: http://gatkforums.broadinstitute.org/gatk/discussion/1601/how-can-i-prepare-a-fasta-file-to-use-as-reference


# Manual


**Dependencies:**
* python 2.7 or higher :     [www.python.org](https://www.python.org/)
* numpy 1.7.0 or higher :    [www.numpy.org](http://www.numpy.org/) 
* scipy 0.14.0 or higher :   [www.scipy.org](http://www.scipy.org/)
* GATK 2.3 or higher :       [www.broadinstitute.org/gatk/download](http://www.broadinstitute.org/gatk/download/)
* java :                     [http://java.com](http://java.com/en/download/)


**Setting environmental variables:**   
To use Conpair you need to set 2 environmental variables and a PYTHONPATH variable (e.g. by adding following lines to your .bashrc file):  
```
export CONPAIR_DIR=/your/path/to/CONPAIR  
export GATK_JAR=/your/path/to/GenomeAnalysisTK.jar

export PYTHONPATH=${PYTHONPATH}:/your/path/to/CONPAIR/modules/
```
**Default reference genome:**

To avoid specifying the reference file every time you run Conpair, please make sure that you have the following files in the specified directory:  
`/your/path/to/CONPAIR/data/genomes/human_g1k_v37.fa`  
`/your/path/to/CONPAIR/data/genomes/human_g1k_v37.fa.fai`  
`/your/path/to/CONPAIR/data/genomes/human_g1k_v37.dict`
<br/>

**Most common usage and additional options:**   
To generate pileups (GATK required):
```

${CONPAIR_DIR}/scripts/run_gatk_pileup_for_sample.py -B TUMOR_bam -O TUMOR_pileup
${CONPAIR_DIR}/scripts/run_gatk_pileup_for_sample.py -B NORMAL_bam -O NORMAL_pileup


Optional:
--help                              show help message and exit
--reference REFERENCE               reference genome in the fasta format, two additional files (.fai, .dict) located in the same directory as the fasta file are required. You may choose to avoid specifying the reference by following the steps in the "default reference genome" section above.
--markers MARKERS                   the set of preselected genomic positions in the BED format. Default: ${CONPAIR_DIR}/data/markers/GRCh37.autosomes.phase3_shapeit2_mvncall_integrated.20130502.SNV.genotype.sselect_v4_MAF_0.4_LD_0.8.bed
--conpair_dir CONPAIR_DIR           path to ${CONPAIR_DIR}
--gatk GATK                         path to GATK JAR [$GATK by default]
--java JAVA                         path to JAVA [java by default]
--temp_dir_java TEMP_DIR_JAVA       java temporary directory to set -Djava.io.tmpdir
--xmx_java  XMX_JAVA                Xmx java memory setting [default: 12g]

```
Verifying concordance between two samples (tumor and normal):
```  

${CONPAIR_DIR}/scripts/verify_concordance.py -T TUMOR_pileup -N NORMAL_pileup


Optional:
--help                              show help message and exit
--outfile OUTFILE                   write output to OUTFILE
--normal_homozygous_markers_only    use only normal homozygous positions to calculate concordance between TUMOR and NORMAL 
--min_cov MIN_COV                   require min of MIN_COV in both TUMOR and NORMAL to use the marker
--min_mapping_quality MIN_MAP_QUAL  do not use reads with mapping qual below MIN_MAP_QUAL [default: 10]
--min_base_quality  MIN_BASE_QUAL   do not use reads with base qual below MIN_BASE_QUAL of a specified position [default: 20]
--markers MARKERS                   the set of preselected genomic positions in the TXT format. Default: ${CONPAIR_DIR}/data/markers/GRCh37.autosomes.phase3_shapeit2_mvncall_integrated.20130502.SNV.genotype.sselect_v4_MAF_0.4_LD_0.8.txt

```
Estimating contamination level in both the tumor and the normal:
```

${CONPAIR_DIR}/scripts/estimate_tumor_normal_contamination.py -T TUMOR_pileup -N NORMAL_pileup


Optional:
--help                              show help message and exit
--outfile OUTFILE                   write output to OUTFILE
--min_mapping_quality MIN_MAP_QUAL  do not use reads with mapping qual below MIN_MAP_QUAL [default: 10] 
--markers MARKERS                   the set of preselected genomic positions in the TXT format. Default: ${CONPAIR_DIR}/data/markers/GRCh37.autosomes.phase3_shapeit2_mvncall_integrated.20130502.SNV.genotype.sselect_v4_MAF_0.4_LD_0.8.txt
--conpair_dir CONPAIR_DIR           path to ${CONPAIR_DIR}
--grid  GRID                        grid interval [default: 0.01]

```
# Output files
**Pileup**  
An example of a pileup file (10 first lines) can be viewed here: ([`pileup.txt`](https://github.com/nygenome/Conpair/blob/master/data/example/pileup/NA12878_normal40x.gatk.pileup.10lines.txt)).

**Concordance**  
An example of a concordance file can be viewed here: ([`concordance.txt`](https://github.com/nygenome/Conpair/blob/master/data/example/concordance/NA12878_tumor80x--NA12878_normal40x.concordance.txt)). 

**Contamination**  
An example of a concordance file can be viewed here: ([`contamination.txt`](https://github.com/nygenome/Conpair/blob/master/data/example/contamination/NA12878_tumor80x--NA12878_normal40x.contamination.txt)). 


# Interpretation  
**Concordance**  (一致性)

To eliminate the effect of copy number variation on the concordance levels, we recommend using the -H flag. 
If two samples are concordant the expected concordance level should be close to 99-100%.  

为了消除拷贝数对一致性水`concordance`平的影响，建议使用 -H 标志。如果两个样本一致，则预期一致性水平应接近 99-100%。

For discordant samples concordance level should be close to 40%.  

对于不一致的样本，一致性水平应接近 **40%**。

You can observe slighly lower concordance (80-99%) in presence of contamination and/or copy number changes (if the -H option wasn't used) in at least one of the samples.   

在存在污染和/或拷贝数变化（如果未使用 -H 选项）的情况下，您可以在至少一个样品中观察到略低的一致性 （**80-99%**）。
**Contamination**(污染)

Even a very low contamination level (such as 0.5%) in the tumor sample will have a severe effect on calling somatic mutations, resulting in decreased specificity. Cross-individual contamination in the normal sample usually has a milder effect on somatic calling.

即使`肿瘤样本`中的**污染水平非常低（例如 0.5%）**也会对calling 体细胞突变产生严重影响，从而导致特异性降低。`正常样本`中的跨个体污染通常对体细胞调用的**影响较轻**。

# 郑

修改参考：[Issues · nygenome/Conpair](https://github.com/nygenome/Conpair/issues/11)

## 特点

- 支持的基因组

GRCh37, GRCh38, GRCm38

- 出发文件的选择

就测试的样本而言，sort.bam和bsqr.bam的结果相差不大。

- Python版本

python2与python3均可适用，包括新增的脚本`run_filter_plieup.py` 

- 消耗时间

整体耗时较短，使用数分钟就可以完成Conpair分析

## 使用注意事项

可以使用环境变量指定`CONPAIR_DIR` ，`GATK_JAR` ，`PYTHONPATH` 来指定其中的选项。这里不使用环境变量的方法去指定这些内容。该方法在迁移时可能会比较麻烦。

### run_gatk_pileup_for_sample.py

- 新增了`--gatk4` 选项，表明调用的是GATK4

- `--conpair_dir`  软件`Conpair` 的路径

- `--reference` 参考基因组，`genome.fa`，`genome.dict`，`genome.fa.fai`

- `--markers` 默认调取GRCh37基因组版本`markers` 

  `GRCh37.autosomes.phase3_shapeit2_mvncall_integrated.20130502.SNV.genotype.sselect_v4_MAF_0.4_LD_0.8.bed`

其它文件在`Conpair` 文件夹下。

- `--gatk` 

对于GATK3，需要指定`GenomeAnalysisTK.jar` ，需要正常指定`--java`

对于GATK4，（方法一）如果使用`conda` 时想指定`*.jar`，需要参考路径`~/miniconda3/envs/rna/share/gatk4-4.5.0.0-0/gatk-package-4.5.0.0-local.jar`，，需要正常指定`--java`；（方法二）如果想指定`gatk`文件，也是没有问题的，参考`~/miniconda3/envs/rna/bin/gatk ` ，此时`--java` 内容不会被调用。

- `--xmx_java` java内存设置，默认为`12g`

### run_filter_pileup.py

- 功能：过滤上一步产生的pileup文件，用于`verify_concordance.py` 的输入文件。
- `--pileup` 指定上一步的`pileup` 文件结果
- `--outfile` 用于`verify_concordance.py` 的结果文件。默认为`<pileup>4concordance`
- `--emptyfile ` 被过滤掉的位点。默认为`<pileup>.empty`，不过滤则会报错。

### verify_concordance.py

- `--tumor_pileup` 和`--normal_pileup` 为过滤后的pileup文件
- `--conpair_dir` 软件`Conpair` 的路径
- `--outfile` 一致性`concordance` 结果文件
- `--markers` 默认调用GRCh37版本基因组markers，`GRCh37.autosomes.phase3_shapeit2_mvncall_integrated.20130502.SNV.genotype.sselect_v4_MAF_0.4_LD_0.8.txt`
- `--min_cov` 最小覆盖度，默认为10
- `--min_mapping_quality` 最小比对质量，默认为10
- `--min_base_quality` 最小碱基质量，默认为20
- `-H,--normal_homozygous_markers_only` 移除拷贝数变对变异的影响。【作者推荐添加，添加后一致性concordance会增高】

#### 不添加`-H`结果

```sh
$ cat WES_FD_1_concordance.txt
Concordance: 71.68%
Based on 6092/7353 markers (coverage per marker threshold : 10 reads)
Minimum mappinq quality: 10
Minimum base quality: 20
```

#### 添加`-H`结果

```sh
$ cat WES_FD_1_concordance_H.txt
Concordance: 96.33%
Based on 3431/7353 markers (coverage per marker threshold : 10 reads)
Minimum mappinq quality: 10
Minimum base quality: 20
```

### estimate_tumor_normal_contamination.py

- `--tumor_pileup` 和`--normal_pileup` 为过滤后的pileup文件
- `--conpair_dir` 软件`Conpair` 的路径
- `--markers` 默认调用GRCh37版本基因组markers，`GRCh37.autosomes.phase3_shapeit2_mvncall_integrated.20130502.SNV.genotype.sselect_v4_MAF_0.4_LD_0.8.txt`
- `--outfile` 污染检测结果文件
- `--grid` 默认值`0.01` 
- `--min_mapping_quality` 最小比对质量，默认10

#### 结果

```sh
$ cat WES_FD_1_contamination.txt
Normal sample contamination level: 0.291%
Tumor sample contamination level: 0.351%
```

## 使用python2【建议】

### run_gatk_pileup_for_sample.py

Tumor：方法一运行run_gatk_pileup_for_sample.py

- --gatk 指定`local.jar` 

这里参考了

```sh
~/miniconda3/envs/py27/bin/python ~/biosoft/Conpair/scripts/run_gatk_pileup_for_sample.py \
--gatk4 \
--bam WES_FD_T_1.bsqr.bam \
--outfile WES_FD_T_1_pileup \
--conpair_dir ~/biosoft/Conpair/ \
--reference ~/db/ref/ucsc-human-hg38/hg38.fa \
--markers ~/biosoft/Conpair/data/markers/GRCh38.autosomes.phase3_shapeit2_mvncall_integrated.20130502.SNV.genotype.sselect_v4_MAF_0.4_LD_0.8.liftover.bed \
--gatk ~/miniconda3/envs/rna/share/gatk4-4.5.0.0-0/gatk-package-4.5.0.0-local.jar \
--java ~/miniconda3/envs/rna/bin/java \
--xmx_java 20g > WES_FD_T_1_gatk_pileup.log 2>&1

~/miniconda3/envs/py27/bin/python ~/biosoft/Conpair/scripts/run_filter_pileup.py --pileup WES_FD_T_1_pileup --outfile WES_FD_T_1_pileup4concordance  --emptyfile WES_FD_T_1_pileup.empty
# awk -F" " '{if($5!="") print $0}'  WES_FD_T_1_pileup >  WES_FD_T_1_pileup4concordance

```

Normal：方法二 运行run_gatk_pileup_for_sample.py

- --gatk 指定最后为`gatk` ，而不是`***.jar`

- `--java` 该方法下java随意写，因为实际运行时该内容并未调用

```sh
~/miniconda3/envs/py27/bin/python ~/biosoft/Conpair/scripts/run_gatk_pileup_for_sample.py \
--gatk4 \
--bam WES_FD_N_1.bsqr.bam \
--outfile WES_FD_N_1_pileup \
--conpair_dir ~/biosoft/Conpair/ \
--reference ~/db/ref/ucsc-human-hg38/hg38.fa \
--markers ~/biosoft/Conpair/data/markers/GRCh38.autosomes.phase3_shapeit2_mvncall_integrated.20130502.SNV.genotype.sselect_v4_MAF_0.4_LD_0.8.liftover.bed \
--gatk ~/miniconda3/envs/rna/bin/gatk \
--java "" \
--xmx_java 20g > WES_FD_N_1_gatk_pileup.log 2>&1

~/miniconda3/envs/py27/bin/python ~/biosoft/Conpair/scripts/run_filter_pileup.py --pileup WES_FD_N_1_pileup --outfile WES_FD_N_1_pileup4concordance  --emptyfile WES_FD_N_1_pileup.empty
# awk -F" " '{if($5!="") print $0}'  WES_FD_N_1_pileup >  WES_FD_N_1_pileup4concordance
```

### verify_concordance.py

```sh
~/miniconda3/envs/py27/bin/python ~/biosoft/Conpair/scripts/verify_concordance.py \
-T WES_FD_T_1_pileup4concordance \
-N WES_FD_N_1_pileup4concordance \
--outfile WES_FD_1_concordance.txt \
--markers ~/biosoft/Conpair/data/markers/GRCh38.autosomes.phase3_shapeit2_mvncall_integrated.20130502.SNV.genotype.sselect_v4_MAF_0.4_LD_0.8.liftover.txt \
--conpair_dir ~/biosoft/Conpair/ \
--min_cov 10 \
--min_mapping_quality 10 \
--min_base_quality 20 > WES_FD_1_verify_concordance.log 2>&1
```

### estimate_tumor_normal_contamination.py

```sh
~/miniconda3/envs/py27/bin/python ~/biosoft/Conpair/scripts/estimate_tumor_normal_contamination.py \
-T WES_FD_T_1_pileup4concordance \
-N WES_FD_N_1_pileup4concordance \
--outfile WES_FD_1_contamination.txt \
--markers ~/biosoft/Conpair/data/markers/GRCh38.autosomes.phase3_shapeit2_mvncall_integrated.20130502.SNV.genotype.sselect_v4_MAF_0.4_LD_0.8.liftover.txt \
--grid 0.01 \
--min_mapping_quality 10 > WES_FD_1_verify_contamination.log 2>&1
# 14:11 开始
# 14:16 结束
```

## 使用python3

### run_gatk_pileup_for_sample.py

Tumor：方法一运行run_gatk_pileup_for_sample.py

````sh
~/miniconda3/envs/rna/bin/python ~/biosoft/Conpair/scripts/run_gatk_pileup_for_sample.py \
--gatk4 \
--bam WES_FD_T_1.bsqr.bam \
--outfile WES_FD_T_1_pileup \
--conpair_dir ~/biosoft/Conpair/ \
--reference ~/db/ref/ucsc-human-hg38/hg38.fa \
--markers ~/biosoft/Conpair/data/markers/GRCh38.autosomes.phase3_shapeit2_mvncall_integrated.20130502.SNV.genotype.sselect_v4_MAF_0.4_LD_0.8.liftover.bed \
--gatk ~/miniconda3/envs/rna/share/gatk4-4.5.0.0-0/gatk-package-4.5.0.0-local.jar \
--java ~/miniconda3/envs/rna/bin/java \
--xmx_java 20g > WES_FD_T_1_gatk_pileup.log 2>&1

~/miniconda3/envs/rna/bin/python ~/biosoft/Conpair/scripts/run_filter_plieup.py --pileup WES_FD_T_1_pileup --outfile WES_FD_T_1_pileup4concordance  --emptyfile WES_FD_T_1_pileup.empty
# awk -F" " '{if($5!="") print $0}'  WES_FD_T_1_pileup >  WES_FD_T_1_pileup4concordance

````

Normal：方法二 运行run_gatk_pileup_for_sample.py

```sh
~/miniconda3/envs/rna/bin/python ~/biosoft/Conpair/scripts/run_gatk_pileup_for_sample.py \
--gatk4 \
--bam WES_FD_N_1.bsqr.bam \
--outfile WES_FD_N_1_pileup \
--conpair_dir ~/biosoft/Conpair/ \
--reference ~/db/ref/ucsc-human-hg38/hg38.fa \
--markers ~/biosoft/Conpair/data/markers/GRCh38.autosomes.phase3_shapeit2_mvncall_integrated.20130502.SNV.genotype.sselect_v4_MAF_0.4_LD_0.8.liftover.bed \
--gatk ~/miniconda3/envs/rna/bin/gatk \
--java "" \
--xmx_java 20g > WES_FD_N_1_gatk_pileup.log 2>&1

~/miniconda3/envs/rna/bin/python ~/biosoft/Conpair/scripts/run_filter_plieup.py --pileup WES_FD_N_1_pileup --outfile WES_FD_N_1_pileup4concordance  --emptyfile WES_FD_N_1_pileup.empty
# awk -F" " '{if($5!="") print $0}'  WES_FD_N_1_pileup >  WES_FD_N_1_pileup4concordance
```

### verify_concordance.py

- 不添加`--normal_homozygous_markers_only`

```sh
~/miniconda3/envs/rna/bin/python ~/biosoft/Conpair/scripts/verify_concordance.py \
-T WES_FD_T_1_pileup4concordance \
-N WES_FD_N_1_pileup4concordance \
--outfile WES_FD_1_concordance.txt \
--markers ~/biosoft/Conpair/data/markers/GRCh38.autosomes.phase3_shapeit2_mvncall_integrated.20130502.SNV.genotype.sselect_v4_MAF_0.4_LD_0.8.liftover.txt \
--conpair_dir ~/biosoft/Conpair/ \
--min_cov 10 \
--min_mapping_quality 10 \
--min_base_quality 20 \
--normal_homozygous_markers_only > WES_FD_1_verify_concordance.log 2>&1
```



- 添加`--normal_homozygous_markers_only`

为了消除拷贝数对一致性水`concordance`平的影响，建议使用 -H 标志。

```sh
~/miniconda3/envs/rna/bin/python ~/biosoft/Conpair/scripts/verify_concordance.py \
-T WES_FD_T_1_pileup4concordance \
-N WES_FD_N_1_pileup4concordance \
--outfile WES_FD_1_concordance_H.txt \
--markers ~/biosoft/Conpair/data/markers/GRCh38.autosomes.phase3_shapeit2_mvncall_integrated.20130502.SNV.genotype.sselect_v4_MAF_0.4_LD_0.8.liftover.txt \
--conpair_dir ~/biosoft/Conpair/ \
--min_cov 10 \
--min_mapping_quality 10 \
--min_base_quality 20 \
--normal_homozygous_markers_only > WES_FD_1_verify_concordance.log 2>&1
```

### estimate_tumor_normal_contamination.py

```sh
~/miniconda3/envs/rna/bin/python ~/biosoft/Conpair/scripts/estimate_tumor_normal_contamination.py \
-T WES_FD_T_1_pileup4concordance \
-N WES_FD_N_1_pileup4concordance \
--outfile WES_FD_1_contamination.txt \
--markers ~/biosoft/Conpair/data/markers/GRCh38.autosomes.phase3_shapeit2_mvncall_integrated.20130502.SNV.genotype.sselect_v4_MAF_0.4_LD_0.8.liftover.txt \
--grid 0.01 \
--min_mapping_quality 10 > WES_FD_1_verify_contamination.log 2>&1
# 16:34 开始
# 16:37 结束
```

## 帮助文档

### run_gatk_pileup_for_sample.py

```sh
$ ~/miniconda3/envs/py27/bin/python ~/biosoft/Conpair/scripts/run_gatk_pileup_for_sample.py
Usage: run_gatk_pileup_for_sample.py [options]

Program to run GATK Pileup on a single sample

Options:
  --version             show program's version number and exit
  -h, --help            show this help message and exit
  -B BAM, --bam=BAM     BAMFILE [mandatory field]
  -O OUTFILE, --outfile=OUTFILE
                        OUTPUT FILE (PILEUP) [mandatory field]
  -D CONPAIR_DIR, --conpair_dir=CONPAIR_DIR
                        CONPAIR DIR [$CONPAIR_DIR by default]
  -R REFERENCE, --reference=REFERENCE
                        REFERENCE GENOME [GRCh37 by default]
  -M MARKERS, --markers=MARKERS
                        MARKER FILE [GRCh37-default]
  -G GATK, --gatk=GATK  GATK JAR [$GATK by default]
  -J JAVA, --java=JAVA  PATH to JAVA [java by default]
  -t TEMP_DIR_JAVA, --temp_dir_java=TEMP_DIR_JAVA
                        temporary directory to set -Djava.io.tmpdir
  -m XMX_JAVA, --xmx_java=XMX_JAVA
                        Xmx java memory setting [default: 12g]
  --remove_chr_prefix   REMOVE CHR PREFIX FROM THE CHROMOSOME COLUMN IN THE
                        OUTPUT FILE [false by default]
  --gatk4               Use GATK4 instead of GATK3 or even lower versions

```

### run_filter_pileup.py

```sh
$ ~/miniconda3/envs/py27/bin/python ~/biosoft/Conpair/scripts/run_filter_pileup.py
Usage: run_filter_pileup.py [options]

Program to filter GATK Pileup on a single sample

Options:
  --version             show program's version number and exit
  -h, --help            show this help message and exit
  -i PILEUP, --pileup=PILEUP
                        pileup file [ required ]
  -o OUTFILE, --outfile=OUTFILE
                        filtered pileup file for script verify_concordance.py
                        [ <pileup>4concordance ]
  -e EMPTYFILE, --emptyfile=EMPTYFILE
                        empty pileup file [ <pileup>.empty ]
```

### verify_concordance.py

```sh
$  ~/miniconda3/envs/py27/bin/python ~/biosoft/Conpair/scripts/verify_concordance.py
Usage: verify_concordance.py [options]

Program to verify tumor-normal sample concordance

Options:
  --version             show program's version number and exit
  -h, --help            show this help message and exit
  -T TUMOR_PILEUP, --tumor_pileup=TUMOR_PILEUP
                        TUMOR PILEUP FILE [mandatory field]
  -N NORMAL_PILEUP, --normal_pileup=NORMAL_PILEUP
                        NORMAL PILEUP FILE [mandatory field]
  -D CONPAIR_DIR, --conpair_dir=CONPAIR_DIR
                        CONPAIR DIR [$CONPAIR_DIR by default]
  -M MARKERS, --markers=MARKERS
                        MARKER FILE [Conpair-GRCh37-default]
  -C MIN_COV, --min_cov=MIN_COV
                        MIN COVERAGE TO CALL GENOTYPE [default: 10]
  -O OUTFILE, --outfile=OUTFILE
                        TXT OUTPUT FILE [stdout by default]
  -Q MIN_MAPPING_QUALITY, --min_mapping_quality=MIN_MAPPING_QUALITY
                        MIN MAPPING QUALITY [default: 10]
  -B MIN_BASE_QUALITY, --min_base_quality=MIN_BASE_QUALITY
                        MIN BASE QUALITY [default: 20]
  -H, --normal_homozygous_markers_only
                        USE ONLY MARKERS THAT ARE HOMOZYGOUS IN THE NORMAL
                        SAMPLE (concordance will not be affected by CNV)
```

### estimate_tumor_normal_contamination.py

```sh
$ ~/miniconda3/envs/py27/bin/python ~/biosoft/Conpair/scripts/estimate_tumor_normal_contamination.py
Usage: estimate_tumor_normal_contamination.py [options]

Program to estimate tumor-normal sample contamination

Options:
  --version             show program's version number and exit
  -h, --help            show this help message and exit
  -T TUMOR_PILEUP, --tumor_pileup=TUMOR_PILEUP
                        TUMOR PILEUP FILE [mandatory field]
  -N NORMAL_PILEUP, --normal_pileup=NORMAL_PILEUP
                        NORMAL PILEUP FILE [mandatory field]
  -D CONPAIR_DIR, --conpair_dir=CONPAIR_DIR
                        CONPAIR DIR [default: $CONPAIR_DIR]
  -M MARKERS, --markers=MARKERS
                        MARKER FILE [default: markers for GRCh37 from
                        $CONPAIR_DIR/data/markers/ ]
  -O OUTFILE, --outfile=OUTFILE
                        TXT OUTPUT FILE [default: stdout]
  -G GRID, --grid=GRID  GRID INTERVAL [default: 0.01]
  -Q MIN_MAPPING_QUALITY, --min_mapping_quality=MIN_MAPPING_QUALITY
                        MIN MAPPING QUALITY [default: 10]
```



## 报错与解决

### 报错1

命令

```sh
~/miniconda3/envs/py27/bin/python ~/biosoft/Conpair/scripts/verify_concordance.py \
-T WES_FD_T_1_pileup \
-N WES_FD_N_1_pileup \
--outfile WES_FD_1_concordance.xls \
--markers ~/biosoft/Conpair/data/markers/GRCh38.autosomes.phase3_shapeit2_mvncall_integrated.20130502.SNV.genotype.sselect_v4_MAF_0.4_LD_0.8.liftover.txt \
--conpair_dir ~/biosoft/Conpair/ \
--min_cov 10 \
--min_mapping_quality 10 \
--min_base_quality 20 > WES_FD_1_verify_concordance.log 2>&1
```

报错

```sh
Traceback (most recent call last):
  File "~/biosoft/Conpair/scripts/verify_concordance.py", line 68, in <module>
    Normal_genotype_likelihoods = genotype_likelihoods_for_markers(Markers, opts.normal_pileup, min_map_quality=MMQ, min_base_quality=MBQ)
  File "~/biosoft/Conpair/modules/ContaminationMarker.py", line 106, in genotype_likelihoods_for_markers
    pileup = parse_mpileup_line(line, min_map_quality=min_map_quality, min_base_quality=min_base_quality)
  File "~/biosoft/Conpair/modules/ContaminationMarker.py", line 71, in parse_mpileup_line
    baseQs = baseQ2int(line[4])
IndexError: list index out of range

```

解决

- 原因：增加打印碱基质量的列表`baseQs`，发现`WES_FD_N_1_pileup` 中第1042行的质量值是空（实际是第4，5列为空），造成列表结果为空而报错。

- 解决思路：对`pileup`质量值列进行检查，去掉该行为空的值。
