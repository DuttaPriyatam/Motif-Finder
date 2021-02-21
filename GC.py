#importing important libraries
import pandas as pd
import numpy as np
from window_slider import Slider

#reading locations from the file
loc_file=open('plus_try.txt', 'r')
lines= loc_file.readlines()
loc_file.close()
pos=[]
for num in lines:
    pos.append(int(num))

#reading the sequences
genome=open('Pfal_try.txt', 'r')

#window size and overlap: parameters for Widnow Slider
size=2
overlap=1

#Declaring arrays for the dictionary
windows=[]
dnts=[]
G_count=[]
C_count=[]
tot_count=[]

#Accessing each position in the location file and seeking the position in our genome file
for pos1 in pos:
    genome.seek(pos1)
    frame=genome.read(1000)
    res=list(frame)
    res_arr=np.array(res)

#Using Window Slider to slide the window
    slider=Slider(size, overlap)
    slider.fit(res_arr)
    i=0

#While the genome does not end/ reach the end the following tasks are performed
    while True:
        data=slider.slide()
        Window=str(pos1+i)+"-"+str(pos1+i+1)
        windows.append(Window)
        dnt=''.join(data)
        dnts.append(dnt)

#Finding the G,C and toal GC for each window
        G_pat='G'
        C_pat='C'

        Gct=dnt.count(G_pat)
        Cct=dnt.count(C_pat)

        G_count.append(Gct)
        C_count.append(Cct)
        tot_count.append(Gct+Cct)
        i=i+1

#If the end of the genome is reached the dictionary is updated for each location sequence.
        if slider.reached_end_of_list():
            dict1= {'Windows':windows,'Dnts':dnts, 'G' : G_count, 'C' : C_count, 'Total':tot_count}
            df=pd.DataFrame.from_dict(dict1)
            df=df.transpose()
            print(df)
            str1="gene"+str(pos1)
            df.to_csv(str1)
            windows=[]
            dnts=[]
            G_count=[]
            C_count=[]
            tot_count=[]
            dict1.clear()
            break
#Clearing the dictionary in the end so a new dictionary can be updated
    dict1.clear()
