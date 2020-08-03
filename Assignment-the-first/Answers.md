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
   
![image](https://github.com/2020-bgmp/demultiplexing-nikkiszc/blob/master/plot_Index1.png) 
![image](https://github.com/2020-bgmp/demultiplexing-nikkiszc/blob/master/plot_Index2.png)
![image](https://github.com/2020-bgmp/demultiplexing-nikkiszc/blob/master/plot_Read1.png) 
![image](https://github.com/2020-bgmp/demultiplexing-nikkiszc/blob/master/plot_Read2.png)
    
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
