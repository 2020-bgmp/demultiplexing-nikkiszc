# Assignment the First

## Part 1
1. Be sure to upload your Python script.

| File name | label |
|---|---|
| 1294_S1_L008_R1_001.fastq.gz | read1 |
| 1294_S1_L008_R2_001.fastq.gz | index1 |
| 1294_S1_L008_R3_001.fastq.gz | index2 |
| 1294_S1_L008_R4_001.fastq.gz | read2 |

2. Per-base NT distribution

   i. Turn in the 4 histograms.
   
![image](https://github.com/2020-bgmp/demultiplexing-nikkiszc/blob/master/plot_Index1.png) 
![image](https://github.com/2020-bgmp/demultiplexing-nikkiszc/blob/master/plot_Index2.png)
![image](https://github.com/2020-bgmp/demultiplexing-nikkiszc/blob/master/plot_Read1.png) 
![image](https://github.com/2020-bgmp/demultiplexing-nikkiszc/blob/master/plot_Read2.png)

   ii. What is a good quality score cutoff for index reads and biological read pairs to utilize 
   for sample identification and downstream analysis, respectively? Justify your answer.

```
After a cursory google search, it seems like a qscore of 30 is a good cutoff. 
This seems reasonable as 30 corresponds to a probability of an error of 0.001, or a 0.1% chance that the base call is incorrect.
This corresponds to an accuracy probability of 99.9%.
In practice, I would test out several options for the qscore cutoff and see what kind of an effect is has on my data.
33 is another contender, as this represents an error probability of 0.0005, or 0.05%, although I don't know if this is too stringent.
```

   iii. How many indexes have undetermined (N) base calls? (Utilize your command line tool knowledge. Submit the command you used. CHALLENGE: use a one-line command)
   
```zcat 1294_S1_L008_R2_001.fastq.gz 1294_S1_L008_R3_001.fastq.gz | sed -n '2~4p' | grep "N" | wc -l```

```7304664```





## Part 2
1. Define the problem

```
We want to de-multiplex the four fastq files we are given, 
meaning that we want to separate different libraries into their own files.
In addition, the forward and reverse reads must be placed in separate files.
Finally, there must be a way to get rid off mismatched or untrustworthy indexes, 
since this obviously affects the outcome of the aligments for each library.
```

2. Describe output

```
The end result should have separate FASTQ files for each of the forward and reverse reads for every library.
Therefore, because there are 24 different libraries, there need to be 48 ouput FASTQ files..
In addition, there need to be a forward and reverse file for both unmatched and unknown reads.
```

3. Upload your [4 input FASTQ files](../TEST-input_FASTQ) and your [4 expected output FASTQ files](../TEST-output_FASTQ).

```See the appropriate folders.```

4. Pseudocode

```See file named pseudocode.txt.```

5. High level functions. For each function, be sure to include:
    1. Description/doc string
    2. Function headers (name and parameters)
    3. Test examples for individual functions
    4. Return statement
    
```This is all included in my pseudocode.```
