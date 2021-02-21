import pandas as pd
import numpy as np
from window_slider import Slider
import getopt
import sys
import time
import re

#Function to take the file names as argument from the command line
def take_arguments():
    argv=sys.argv[1:]
    InputFile=None
    LocationFile=None
    motif=None

    try:
        opts, args=getopt.getopt(argv, "hi:l:m:", ["help", "genome_input_file =", "location_file =", "Motif ="])

    except:
        print("Error in the input")
    allowed_bases=['A', 'T', 'G', 'C']
    for opt, arg in opts:
        if opt in ['-h', '--help']:
            sys.exit()
        elif opt in ['-i', '--genome_input_file']:
            InputFile=arg
        elif opt in ['-l', '--location_file']:
            LocationFile=arg
        elif opt in ['-m', '--motif']:
            motif=arg


    file_list=[InputFile, LocationFile, motif]
    return(file_list)

#Obtain locations from the location file
def get_locations(loc_file):
    loc_file=open(loc_file, 'r')
    lines= loc_file.readlines()
    loc_file.close()
    pos=[]
    for num in lines:
        pos.append(int(num))
    return(pos)

#Function to find particular motif, create file for each gene 
def finding_motifs(gene_positions, gen_file, motif):
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
    print(List)

    for x in List:
        motif_rev="".join(complement.get(base, base) for base in reversed(x))
        reverse_List.append(motif_rev)
    print(reverse_List)
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

    #Finding the start positions, creating 1000 bp frame, and storing it as an array.
    final=0
    for pos in gene_positions:
        genome.seek(pos-1)
        frame=genome.read(1000)
        res=list(frame)
        res_arr=np.array(res)

    #Creating window slider of required size and overlap.
        slider=Slider(size, overlap)
        slider.fit(res_arr)
        i=0

        template=0
        non_template=0
        dict={}

        while True:
            data=slider.slide()
            Window=str(pos+i)+"-"+str(pos+i+(size-1))
            windows.append(Window)
            string=''.join(data)
            strings.append(string)

    #Finding PAM motif and printing frequency.
            count=List.count(string)
            PAM_count.append(count)
            non_template=non_template+count
            
            count_rev=reverse_List.count(string)
            PAM_rev_count.append(count_rev)
            template=template+count_rev
            i=i+1

            if slider.reached_end_of_list():
                dict1= {'Windows':windows,'PAM':strings, motif : PAM_count, reverse : PAM_rev_count}
                df=pd.DataFrame.from_dict(dict1)
                df=df.transpose()
                print(df)
                str1="gene"+str(pos)
                df.to_csv(str1)
                windows=[]
                strings=[]
                PAM_count=[]
                PAM_rev_count=[]
                dict1.clear()
                break
        dict1.clear()
        total_temp=total_temp+template
        total_non_temp=total_non_temp+non_template
    
    count_list=[total_temp, total_non_temp]
    return(count_list)


begin=time.time()
file_list=take_arguments()
gene_positions=get_locations(file_list[1])
count_list=finding_motifs(gene_positions, file_list[0], file_list[2])

total=count_list[0]+count_list[1]
average=total/2901

print("Total PAMS on non-template strand :"+str(count_list[0]))
print("Total PAMS on template strand :"+str(count_list[1]))
print("Total count of PAMs: "+str(total))
print("Average occurrence of PAM :"+str(average))
time.sleep(1)
end=time.time()
print(f"Total runtime of the program is {end - begin}")
