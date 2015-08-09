import sys
import os
import itertools

diseases={}

inputfile=sys.argv[2]
#print inputfile
#inputfile='sample.txt'
outputfile=inputfile.split('/')[-1].replace(".txt","_inference.txt")
#print outputfile
#outputfile="$LIB/sample_input_inference.txt"
#print outputfile
def calculate_prob(finding_list,diseases,dis):
    #print finding_list,diseases
        #print "dis keys=",diseases.keys()
    #for d,dis in enumerate(diseases.keys()):
        numerator=1.0*diseases[dis]["p(d)"]
        denominator=1.0*(1-diseases[dis]["p(d)"])
        #print "Finding probability of disease=",dis
        #output[dis]=0.0
        for i,fin in enumerate(finding_list):
            if(fin=="T"):
                numerator*=diseases[dis]["p(f|d)"][i]
                denominator*=diseases[dis]["p(f|~d)"][i]
                continue
            if(fin=="F"):
                numerator*=(1-diseases[dis]["p(f|d)"][i])
                denominator*=(1-diseases[dis]["p(f|~d)"][i])
                continue
            if(fin=="U"):
                #Ignoring Unknown findings
                continue
        #print "numerator=",numerator
        #print "denominator=",denominator
        if(numerator+denominator)==0:
                prob=float('inf')
        else:
                prob=numerator/(numerator+denominator)
        #print "Prob of dis {0} is {1}".format(dis,prob)
        #output[dis]=round(prob,4)
        #numerator=1.0
        #denominator=1.0
        return round(prob,4)    


def calculate_minmax(dis_prob,finding_list,diseases,output_mnmx,dis):
        min_prob = max_prob = dis_prob
        ucount=finding_list.count("U")
        copy_finding_list=list(finding_list)
        #print "Original finding_list=",finding_list
        #print "Ucount=",ucount
        if(ucount==0):
                output_mnmx[dis]=[]
                output_mnmx[dis].append(str("{:.4f}".format(min_prob)))
                output_mnmx[dis].append(str("{:.4f}".format(max_prob)))
                return
        perm_list=[]
        for i in list(itertools.product(['T','F'],repeat=ucount)):
                perm_list.append(list(i))
        #print perm_list
        
        for i,per in enumerate(perm_list):
                #print "per=",per
                u_no=0
                for j,ele in enumerate(finding_list):
                        if(ele=="U"):
                              copy_finding_list[j]=perm_list[i][u_no]
                              u_no+=1
                #print "Copy list is",copy_finding_list
                temp_prob=calculate_prob(copy_finding_list,diseases,dis)
                if temp_prob<min_prob:
                        min_prob=temp_prob
                if temp_prob>max_prob:
                        max_prob=temp_prob
        output_mnmx[dis]=[]
        output_mnmx[dis].append(str("{:.4f}".format(min_prob)))
        output_mnmx[dis].append(str("{:.4f}".format(max_prob)))

def calculate_incrdecr(dis_prob,finding_list,diseases,output_id,dis):
        min_prob = max_prob = dis_prob
        #print "dis=",dis
        #print "dis_prob=",dis_prob
        ucount=finding_list.count("U")
        copy_finding_list=list(finding_list)
        #print "Original finding_list=",finding_list
        #print "Ucount=",ucount
        if(ucount==0):
                output_id[dis]=[]
                output_id[dis].append('none')
                output_id[dis].append('N')
                output_id[dis].append('none')
                output_id[dis].append('N')
                return
        
        min_index=max_index=-1
        #Check if prob increases or decreases when the test has been done/not done
        for j,ele in enumerate(finding_list):
                copy_finding_list=list(finding_list)
                if(ele=="U"):
                        copy_finding_list[j]='T'
                        
                        #print "Copy list is",copy_finding_list
                        temp_prob=calculate_prob(copy_finding_list,diseases,dis)
                        #print "Temp_prob=",temp_prob
                        if temp_prob<min_prob:
                                min_prob=temp_prob
                                min_index=j
                                min_per='T'
                        #Alphabetical ordering of symptoms if same prob
                        if temp_prob==min_prob and min_index!=-1:
                                x=diseases[dis]["findings"][min_index]
                                y=diseases[dis]["findings"][j]
                                if(x>y):
                                        min_index=j
                        if temp_prob>max_prob:
                                max_prob=temp_prob
                                max_index=j
                                max_per='T'
                        if temp_prob==max_prob and max_index!=-1:
                                x=diseases[dis]["findings"][max_index]
                                y=diseases[dis]["findings"][j]
                                if(x>y):
                                        max_index=j
                        
                        copy_finding_list[j]='F'
                        
                        #print "Copy list is",copy_finding_list
                        temp_prob=calculate_prob(copy_finding_list,diseases,dis)
                        #print "Temp_prob=",temp_prob
                        if temp_prob<min_prob:
                                min_prob=temp_prob
                                min_index=j
                                min_per='F'
                        if temp_prob==min_prob and min_index!=-1:
                                x=diseases[dis]["findings"][min_index]
                                y=diseases[dis]["findings"][j]
                                if(x>y):
                                        min_index=j
                        if temp_prob>max_prob:
                                max_prob=temp_prob
                                max_index=j
                                max_per='F'
                        if temp_prob==max_prob and max_index!=-1:
                                x=diseases[dis]["findings"][max_index]
                                y=diseases[dis]["findings"][j]
                                if(x>y):
                                        max_index=j
                        
        output_id[dis]=[]
        if(max_index!=-1):
                output_id[dis].append(diseases[dis]["findings"][max_index])
                output_id[dis].append(max_per)
        else:
                output_id[dis].append('none')
                output_id[dis].append('N')
        if(min_index!=-1):
                output_id[dis].append(diseases[dis]["findings"][min_index])
                output_id[dis].append(min_per)
        else:        
                output_id[dis].append('none')
                output_id[dis].append('N')
        




        
with open(inputfile, 'r') as f:
        
        #List to preserve the disease correct order to calc. prob. of each Pt.
        disease_list=[]
        line=f.readline()
        no_of_diseases=int(line.split()[0])
        no_of_patients=int(line.split()[1])
        #print "no. of diseases = ",no_of_diseases
        #print "no. of patients = ",no_of_patients
        #Create a dictionary of dictionaries to store all the info regarding disease
        for i in range(1,no_of_diseases+1):
            name,no_of_sym,dis_prob=f.readline().split()
            diseases[name] = {}
            disease_list.append(name)
            diseases[name]["nos"] = int(no_of_sym)
            diseases[name]["p(d)"] = float(dis_prob)
            diseases[name]["findings"]=eval(f.readline())
            diseases[name]["p(f|d)"]=eval(f.readline())
            diseases[name]["p(f|~d)"]=eval(f.readline())
        #print "Dictionary of Diseases= ",diseases
        #print outputfile
        op=open(str(outputfile),'w')
        #For each patient compute probability for 3 questions
        for i in range(1,no_of_patients+1):
            output={}
            output_mnmx={}
            output_id={}
            
            for dis in disease_list:
                finding_list=eval(f.readline())
                #print finding_list
                try:
                	dis_prob=calculate_prob(finding_list,diseases,dis)
                	output[dis]=str("{:.4f}".format(dis_prob))
                	calculate_minmax(dis_prob,finding_list,diseases,output_mnmx,dis)
                	calculate_incrdecr(dis_prob,finding_list,diseases,output_id,dis)
                except:
                    pass
                  
            #print "Patient-{0}:".format(i)
            #print output
            #print output_mnmx
            #print output_id
            
            patient="Patient-"+str(i)+":"  	
            op.write(patient)
            op.write("\n")
            op.write(str(output))
            op.write("\n")
            op.write(str(output_mnmx))
            op.write("\n")
            op.write(str(output_id))
            op.write("\n")
        op.close()   	
  
                



    
