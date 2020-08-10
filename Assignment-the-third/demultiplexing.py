#!/usr/bin/env python

import argparse
from itertools import product
import gzip

def get_args():
    parser = argparse.ArgumentParser(description='Make Qscore distribution graph')
    parser.add_argument("-r1", "--read1_file", type=str, help='Type the input file for read 1')
    parser.add_argument("-r2", "--read2_file", type=str, help='Type the input file for read 2')
    parser.add_argument("-r3", "--read3_file", type=str, help='Type the input file for read 3')
    parser.add_argument("-r4", "--read4_file", type=str, help='Type the input file for read 4')
    parser.add_argument("-i", "--indexes_file", type=str, help='Type the input file for the indexes and index names')
    return parser.parse_args()

parse_args = get_args()
read1_file = parse_args.read1_file
read2_file = parse_args.read2_file
read3_file = parse_args.read3_file
read4_file = parse_args.read4_file
indexes_file = parse_args.indexes_file

#Create dictionary of index names and index sequences
def Create_Index_Name_Dict(fh):
    with open(fh, "r") as inf:
        index_dict = {}
        line_cnt=0
        for line in inf:
            line = line.strip("\n")
            line = line.split("\t")
            line_cnt+=1
            if line_cnt > 1:
                name, index = line[3], line[4]
                index_dict[name]=index
    return index_dict

index_dict = Create_Index_Name_Dict(indexes_file)

#This isn't exactly "permutations" because it also includes the correctly matched pairs of indexes,
#as well as all possible pairs of two. (That's what the repeat=2 means.)
def Create_Permutations_Dictionary(all_indexes):
    Perms_Dict = {}
    for pair in product(all_indexes, repeat = 2):
        Perms_Dict[pair]=0
    return Perms_Dict

Perms_Dict = Create_Permutations_Dictionary(index_dict.keys())

def Reverse_Complement_Function(sequence):
    Rev = ''
    for base in sequence:
        if base == 'A':
            Rev+='T'
        elif base == 'T':
            Rev+='A'
        elif base == 'G':
            Rev+='C'
        elif base == 'C':
            Rev+='G'
        else:
            Rev+='N'
    RevComp = Rev[::-1]
    return RevComp

def Average_Qscore_Function(qscore_line):
    total_score = 0
    qscore_length = 0
    for qscore in qscore_line:
        total_score += ord(qscore)-33
        qscore_length += 1
    average_qscore = total_score/qscore_length
    return average_qscore

def Create_Index_File_Dict(index_dict):
    index_file_dict = {}
    for name in index_dict:
        seq = index_dict[name]
        index_file_dict[seq] = (str(name)+"_fwd.fastq",str(name)+"_rev.fastq")
    return index_file_dict

index_file_dict = Create_Index_File_Dict(index_dict)
#print(index_file_dict)


# Up Next: The big ol' master function.
# This function will split up records from the four Input Files into three new output files: 

# Matched: indexes match, i.e., R2 is the same as the reverse comp of R3
# Unmatched: R2 is NOT the same as the RC of R3
# Unknown: There are Ns in either R2 or R3 or they are not part of the library, or the average Qscore of either index is below 30.

# Because each record in each of the four input files correlate to the same DNA molecule
# (according to line count in the file), read in files simultaneously using zip.

def Categorize_Records_Function(R1,R2,R3,R4):

    Unknown_index_fwd = open('Unknown_index_fwd.fastq', 'a')
    Unknown_index_rev = open('Unknown_index_rev.fastq', 'a')
    Index_hopped_fwd = open('Index_hopped_fwd.fastq', 'a')
    Index_hopped_rev = open('Index_hopped_rev.fastq', 'a')

    #Initialize a list called Current_Record_List that has two lists inside of it,
    #where the first and second lists correspond to the entire record of read1 and read2, respectively.)
    Current_Record_List = [[0,0,0,0],[0,0,0,0]]

    #Simultaneously open up the four input fasta files.
    with gzip.open(R1, "rt") as r1, gzip.open(R2, "rt") as r2, gzip.open(R3, "rt") as r3, gzip.open(R4, "rt") as r4:
        Line_Count = 0
        unknown_indexes_counter = 0
        index_hopped_counter = 0
        matched_indexes_counter = 0
        for r1_line, r2_line, r3_line, r4_line in zip(r1,r2,r3,r4):
            Line_Count += 1
            #I'm keeping the newlines on (for now, at least) since the lines are going to go into a new fastq file.
            #The following if statements store each line of the R1 and R4 record into the Current_Record_List.
            if Line_Count%4 == 1:
                header1 = r1_line.strip()
                header2 = r4_line.strip()
            if Line_Count%4 == 2:
                Current_Record_List[0][1] = r1_line 
                Current_Record_List[1][1] = r4_line 
                r2_seq = r2_line
                r3_seq = r3_line
                r3_rc_seq = Reverse_Complement_Function(r3_seq.strip())
                Current_Record_List[0][0] = (header1 + ' ' + r2_seq.strip() + '-' + r3_rc_seq + '\n')
                Current_Record_List[1][0] = (header2 + ' ' + r2_seq.strip() + '-' + r3_rc_seq + '\n')
                r2_avg_qscore = Average_Qscore_Function(r2_seq.strip())
                r3_avg_qscore = Average_Qscore_Function(r3_seq.strip())
            if Line_Count%4 == 3: 
                Current_Record_List[0][2] = r1_line 
                Current_Record_List[1][2] = r4_line
            if Line_Count%4 == 0: 
                Current_Record_List[0][3] = r1_line 
                Current_Record_List[1][3] = r4_line 
                r2_qscore = r2_line
                r3_qscore = r3_line

                #Now, the Current_Record_List is all filled w the full Read1 and Read2 records.
                #The record now must be moved to the correct file.
                if r2_seq.strip() not in index_dict.values() or r2_avg_qscore < 30 or r3_avg_qscore < 30 or r3_rc_seq not in index_dict.values():
                    #The record is stored in a list, so the lines of each 4-line record must be written one by one.
                    for x in range(4):
                        Unknown_index_fwd.write(str(Current_Record_List[0][x]))
                        Unknown_index_rev.write(str(Current_Record_List[1][x]))
                    unknown_indexes_counter +=1
                
                #If both indexes are valid, but Index1 and Index2 RC don't match (i.e. index hopping):
                elif r2_seq.strip() != r3_rc_seq:
                    index_hopped_counter+=1
                    for x in range(4):
                        Index_hopped_fwd.write(str(Current_Record_List[0][x]))
                        Index_hopped_rev.write(str(Current_Record_List[1][x]))
                    #Adding to Permutations Dictionary
                    index1 = list(index_dict.keys())[list(index_dict.values()).index(r2_seq.strip())]
                    index2 = list(index_dict.keys())[list(index_dict.values()).index(r3_rc_seq)]
                    pair = (index1,index2)
                    Perms_Dict[pair]+=1

                #The only other option now is for the two indexes to be from the index library, and they match.
                else:
                    matched_indexes_counter+=1
                    files = index_file_dict[r2_seq.strip()]
                    fwd_file = open(files[0], 'a')
                    rev_file = open(files[1], 'a')
                    for x in range(4):
                        fwd_file.write(str(Current_Record_List[0][x]))
                        rev_file.write(str(Current_Record_List[1][x]))
                    #Adding to Permutations Dictionary
                    #The following line works like this: key_list[value_list.index(value)]
                    #It lets you get the key from a value in a dictionary, since I want the permutations dictionary 
                    index1 = list(index_dict.keys())[list(index_dict.values()).index(r2_seq.strip())]
                    pair = (index1, index1)
                    Perms_Dict[pair]+=1

    #Close all the files that I opened.
    for fwd_rev in index_file_dict.values():
        fwd = open(fwd_rev[0], "a")
        rev = open(fwd_rev[1], "a")
        fwd.close()
        rev.close()
    Unknown_index_fwd.close()
    Unknown_index_rev.close()
    Index_hopped_fwd.close()
    Index_hopped_rev.close()

    return unknown_indexes_counter, index_hopped_counter, matched_indexes_counter, Perms_Dict

unknown_indexes_counter, index_hopped_counter, matched_indexes_counter, Perms_Dict = Categorize_Records_Function(read1_file,read2_file,read3_file,read4_file)

def Create_Output_File():
    total_reads = unknown_indexes_counter + index_hopped_counter + matched_indexes_counter

    with open("Output_file.txt", "w") as output_fh:
        output_fh.write("\nThis is the output file for the script demultiplexing.py.\n")
        output_fh.write("Number of Matched Index Pairs:" + "\t" + str(matched_indexes_counter) + "\t" + "(" + str(matched_indexes_counter/total_reads*100) + "%)\n")
        output_fh.write("Number of Index-Hopped Index Pairs:" + "\t" + str(index_hopped_counter) + "\t" + "(" + str(index_hopped_counter/total_reads*100) + "%)\n")
        output_fh.write("Number of Unknown/Low Quality Index Pairs:" + "\t" + str(unknown_indexes_counter) + "\t" + "(" + str(unknown_indexes_counter/total_reads*100) + "%)\n")
        output_fh.write("Number of Total Index Pairs:" + "\t" + str(total_reads) + "\n")
        output_fh.write("\nCounters for Matched Index Pairs\n\nIndex Pair\tCount\tPercentage out of Matched Reads\n")

        for key, value in sorted(Perms_Dict.items()):
            if key[0] == key[1]:
                output_fh.write(str(key)+"\t"+str(value) + "\t" + str(value/matched_indexes_counter*100) + "%\n")

        output_fh.write("\nCounters for All Index Pairs\n\nIndex Pair\tCount\tPercentage out of All Reads\n")

        for key, value in sorted(Perms_Dict.items()):
            output_fh.write(str(key)+"\t"+str(value) + "\t" + str(value/total_reads*100) + "%\n")

    return

Create_Output_File()
