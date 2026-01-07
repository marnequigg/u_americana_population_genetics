# American Elm Genetic Analysis
Polyploid complexes, like American elm, are fairly common, but there are very few analysis pipelines designed to handle this. Hopefully, this technique will be helpful for other organisms that have similar polyploid problems. So, I want to start off with some basic background on American elm and what we know about the polyploid problem.

We know that American elm can have four different ploidy levels: 2n, 3n, 4n, and 5n. Tetraploids are the most common, followed by the less common diploids. Then triploids and pentaploids are rare, with only one or two of each (supposedly the triploid 'Jefferson' is fertile and very resistant to Dutch elm disease). Diploids and tetraploids are found across the southern part of the range, and even occur in the same stands. There are no published studies attempting to cross the cytotypes, so it seems improbable, but since the triploids exist, it might be possible. We also don't know if the polyploids are autopolyploids (the result of genome duplication) or allopolyploids (the result of hybridization with two distinct subgenomes). After concluding the analysis for this project, we are fairly certain they are allopolyploids, but still don't have confirmation of who the progenitors are.

I am a strong believer in reporting failures as well as successes, so I am going to try to report everything I tried and didn't work. My goal is not to put down any of the programs that didn't work, just let others know what did and did not work well for my dataset. Feel free to reach out to me if you want to chat more about this.

### Is this pipeline applicable to your data?
So, if you find yourself with a polyploid complex and target enrichment data, this pipeline might be useful for you. Also, if you have an allopolyploid with two *very distinct* (so distinct that SNP callers throw away reads from the other sub genomes), then this might be useful for you (check out part 3!). Also relevant, I used a custom amplicon panel, so I have target enrichment data.

## Part 1) Split Cytotypes
A lot of studies just analyze all the different cytotypes together, but this is *statistically illegal* (no judgment here, I almost did this). The problem arises during SNP-calling because the programs available can only handle one cytotype at a time. Even so, a lot of SNP calling programs really struggle with polyploids (as discussed in part 3). I opted to separate the cytotypes and analyze them completely separate. This works well for my dataset since I have two dominant cytotypes (2n and 4n).

To do this, I used a PCA. In order to use a PCA, I had individuals with known ploidy, as confirmed with flow cytometry. The trees with known ploidy are living at the National Arboretum in Maryland. I also tried using k-mer plots, flow cytometry (dried leaves **do not** work), ploidyNGS, nQuire, nQuack and PCA. nQuack was the best, and I think the problem was with how I was implementing it, but a PCA was the most accurate for my data, and the quickest.

The scripts I used for this section are in a folder titled "part1_scripts".
![graphic of pipeline steps](images/2.png)

### Step 1b) Quality Control
Of course, once I got sequence data back, I ran FastQC and MultiQC. This is just to make sure that sequencing was successful. The associated scripts for this are 1b_fastqc.sh and 1b_multiqc.sh.
### Step 2b) Trimming
This step removes the Illumina adapters used for sequencing. The associates scripts are 2b_trimmomatic.sh
### Step 3b) Alignment
The associated script here is 3b_alignment.sh. This starts with adding the programs to the conda environment and indexing the genome. Then, 

## Part 2) Analyze Diploids

## Part 3) Analyze Polyploids
![graphic of pipeline steps](images/3.png)
