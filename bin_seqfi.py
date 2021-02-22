import pandas as pd
import numpy as np
from window_slider import Slider
import getopt
import math
import sys
import time
import re
import os

#Function to take the file names as argument from the command line
def take_arguments():
    argv=sys.argv[1:]
    input_file=None
    gff_file=None
    motif=None
    bins=None

    try:
        opts, args=getopt.getopt(argv, "hi:g:m:b:", ["help", "genome_input_file =", "gff_file =", "Motif =", "Bins ="])

    except:
        print("Error in the input")
    allowed_bases=['A', 'T', 'G', 'C']
    for opt, arg in opts:
        if opt in ['-h', '--help']:
            sys.exit()
        elif opt in ['-i', '--genome_input_file']:
            input_file=arg
        elif opt in ['-g', '--gff_file']:
            gff_file=arg
        elif opt in ['-m', '--motif']:
            motif=arg
        elif opt in ['-b', '--bins']:
            bins=arg


    file_list=[input_file, gff_file, motif, bins]
    return(file_list)

#Obtain start locations from the gff file
def start_get_locations(loc_file):

    cmd='cut -f1 '+loc_file+' > plus_try_str.txt'
    #cmd='cut -f4 '+loc_file+' > plus_try_str.txt'
    os.system(cmd)

    sloc_file=open('plus_try_str.txt', 'r')
    lines= sloc_file.readlines()
    sloc_file.close()

    start_pos=[]
    for num in lines:
        start_pos.append(int(num))
    return(start_pos)

#Obtain end locations from the gff file
def end_get_locations(loc_file):

    cmd='cut -f2 '+loc_file+' > plus_try_ov.txt'
    #cmd='cut -f5 '+loc_file+' > plus_try_ov.txt'
    os.system(cmd)

    eloc_file=open('plus_try_ov.txt', 'r')
    lines= eloc_file.readlines()
    eloc_file.close()

    end_pos=[]
    for num in lines:
        end_pos.append(int(num))
    return(end_pos)

#Obtain gene ids from the gff
def get_ids(loc_file):

    cmd= 'cut -f3 '+loc_file+' | sed -e "s/ID=//g"  > plus_try_genes.txt'
    #cmd= 'cut -f9 '+loc_file+' | sed -e "s/ID=//g" | sed -e "s/;.*//g" > plus_try_genes.txt'
    os.system(cmd)

    gloc_file=open('plus_try_genes.txt', 'r')
    lines= gloc_file.readlines()
    gloc_file.close()

    genes=[]
    for gene in lines:
        genes.append(gene)
    return(genes)

def get_new_end(start_pos, end_pos, nbin):
    eff_new_end=[]
    for start, end in zip(start_pos, end_pos):
        diff=end-start
        rem=diff%int(nbin)
        efgensi=end-rem
        eff_new_end.append(efgensi)

    return(eff_new_end)

#Making bed annotation for the genes, including the bin size for each gene
def make_gene_data(start_pos, end_pos, genes, eff_end, nbins):

    size=[]
    per_bin=[]
    removed=[]
    for start, end, ends in zip(start_pos, end_pos, eff_end):
        diff=ends-start
        size.append(diff)

        rem=end-ends
        removed.append(rem)

        bin_size=diff/20
        per_bin.append(math.ceil(bin_size))

    gene_bin_dict={'Gene':genes, 'Start':start_pos, 'End':end_pos, 'Effective New End': eff_end, 'Bases Removed':removed, 'Size':size, 'Size per Bin':per_bin}

    bin_df=pd.DataFrame.from_dict(gene_bin_dict)
    print(bin_df)
    bin_df.to_csv("Gene_Bin_Data")

    return(per_bin)

#Function to find particular motif, create file for each gene 
def finding_motifs(start_pos, end_pos, gen_id, per_bin, gen_file, motif):
    #Reading the genome sequence file.
    complement = {'A' : 'T', 'C' : 'G', 'G' : 'C', 'T' : 'A', 'R':'Y', 'Y':'R'}
    genome=open(gen_file, 'r')
    
    
    R=['A','G']
    Y=['C','T']
    N=['A','T','G','C']
    List=[]
    reverse_List=[]
    
    #Using regular expression to check for the presence a character and replacing it with allowed nucleotides in its place and making a list of it
    if re.search('N',motif):
        for a in N:
            motif_replaced=motif.replace('N', a)
            List.append(motif_replaced)
    elif re.search('R', motif):
        for a in R:
            motif_replaced=motif.replace('R', a)
            List.append(motif_replaced)
    elif re.search('Y', motif):
        for a in Y:
            motif_replaced=motif.replace('Y', a)
            List.append(motif_replaced)
    else:
        List.append(motif)

    for x in List:
        motif_rev="".join(complement.get(base, base) for base in reversed(x))
        reverse_List.append(motif_rev)
    reverse="".join(complement.get(base, base) for base in reversed(motif))

    size=len(motif)
    overlap=size-1
    total_temp=0
    total_non_temp=0

    #Defining arrays
    windows=[]
    strings=[]
    PAM_count=[]
    PAM_rev_count=[]
    total_count=[]
    average_count=[]
    #Iterating through our genome file with the start position and end position for each gene and using bin size of that particular gene as the iterator value
    for start, end, bin_size, gene in zip(start_pos, end_pos, per_bin, gen_id):
        ct=1
        x=start
        #for x in range(start, end, bin_size):
        while x <= end:
            if ct==1:
                genome.seek(x)
                frame=genome.read(bin_size+2)
            else:
                x=x-2
                genome.seek(x)
                frame=genome.read(bin_size+2)

            res=list(frame)
            res_arr=np.array(res)

        #Creating window slider of required size and overlap.
            slider=Slider(size, overlap)
            slider.fit(res_arr)
            i=0

    #        template=0
     #       non_template=0

            while True:
                data=slider.slide()
                Window=str(x+i)+"-"+str(x+i+(size-1))
                windows.append(Window)
                string=''.join(data)
                strings.append(string)

        #Finding PAM motif and printing frequency.
                count=List.count(string)
                PAM_count.append(count)
                non_template=non_template+count
            
                count_rev=reverse_List.count(string)
                PAM_rev_count.append(count_rev)
   #             template=template+count_rev
                i=i+1

                if slider.reached_end_of_list():
                    dict1= {'Windows':windows,'PAM':strings, motif : PAM_count, reverse : PAM_rev_count}
                    df=pd.DataFrame.from_dict(dict1)
                    df=df.transpose()
                    print(df)
                    str1=str(gene)+"_bin"+str(ct)
                    ct=ct+1
                    df.to_csv(str1)
                    windows=[]
                    strings=[]
                    PAM_count=[]
                    PAM_rev_count=[]
                    dict1.clear()
                    x=x+bin_size+2
                    break
  #          total_temp=total_temp+template
  #          total_non_temp=total_non_temp+non_template
    
    count_list=[total_temp, total_non_temp]
    return(count_list)


begin=time.time()
file_list=take_arguments()
print(file_list)
start_positions=start_get_locations(file_list[1])
print(start_positions)
end_positions=end_get_locations(file_list[1])
print(end_positions)
ids=get_ids(file_list[1])
print(ids)
new_end=get_new_end(start_positions, end_positions, file_list[3])
print(new_end)
bin_size=make_gene_data(start_positions, end_positions, ids, new_end, file_list[3])
print(bin_size)
count_list=finding_motifs(start_positions, new_end, ids, bin_size, file_list[0], file_list[2])

#total=count_list[0]+count_list[1]
#average=total/2901

#print("Total PAMS on non-template strand :"+str(count_list[0]))
#print("Total PAMS on template strand :"+str(count_list[1]))
#print("Total count of PAMs: "+str(total))
#print("Average occurrence of PAM :"+str(average))

#Removing the temporary files that were generated for the job
os.system('rm plus*')

time.sleep(1)
end=time.time()
print(f"Total runtime of the program is {end - begin}")
