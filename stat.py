import numpy as np
import argparse

parser = argparse.ArgumentParser(description = """Obtain average counts for proplyds, either for the bowshock or the surrounding nebula """)
parser.add_argument("--region",type = str,choices=("shell","nebula"),default = "shell",help = "type of region to measure average counts")
filtro = ["Ha","NII","OIII","con"]
proplyd = ["LV2","LV4","LV5","HST1","167-328"]
#emission = ["","_neb"]

cmd_args = parser.parse_args()
region = cmd_args.region

if region == 'shell':
    for fil in filtro[:-1]:
        for p in proplyd:    
            f=open(fil+p+".pix","r")
            text = f.readlines()
            seltext= text[2:] # The first two lines does not provide useful information
	    sample = []
            for line in seltext:
                linerow= line.split('|')[1] # only the information after the "|" is useful
                row = linerow.split()#changing a string into a list
                sample.append(row) # organize everything in a single list
            Sample = np.array(sample) #changing to numpy array
            S = Sample.astype(np.float) #converting elements from strings to floats
            print fil,p,region,np.mean(S),np.std(S) #shows the mean and the standard deviation of sample

else:
    for fil in filtro:
        for p in proplyd:
            for i in range(0,3):
                f=open(fil+p+"_neb-"+str(i)+".pix","r")
                text = f.readlines()
                seltext= text[2:] # The first two lines does not provide useful information
		sample = []
                for line in seltext:
                       linerow= line.split('|')[1] # only the information after the "|" is useful
                       row = linerow.split()#changing a string into a list
                       sample.append(row) # organize everything in a single list
            Sample = np.array(sample) #changing to numpy array
            S = Sample.astype(np.float) #converting elements from strings to floats
            print fil,p,region,np.mean(S),np.std(S) #shows the mean and the standard deviation of sample



