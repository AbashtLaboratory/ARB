#!/usr/bin/env python3.6

######################################Analysis of RNA-Seq Data Behavior Tool (ARB)##############################################
################################################################################################################################
################################################################################################################################

##############################Code was developed in the Abasht Laboratory at the University of Delaware#########################
####################under the supervision of Dr. Behnam Abasht website: http://canr.udel.edu/faculty/behnam-abasht/#############

#Please contact Dr. Behnam Abasht (abasht@udel.edu) with any questions or concerns.

#Description of Program
#ARB is a python based tool for assessment of data behavior and prioritization of differentiallly expressed genes which uses output files from
#Cuffdiff (Trapnell 2013).  


#Basic Setup of Program

#Part 1. Parsing FPKM Data
#The program first parses through the read_groups.info and genes.read_group_tracking file 
#produced by Cuffdiff to create a more standardized file of FPKM values. The parsed file is 
#then filtered for statistically significant genes using information found in Cuffdiff’s gene_exp.diff file

#Part 2. Inter-Group Separation Index (IGSI) Analysis
#The parsed FPKM data is then used to analyze each gene’s overall data behavior among samples between group 1 
#and group 2 (also known as conditions in Cuffdiff) by calculating the inter-group separation index (IGSI). To 
#calculate IGSI, the program determines if the lowest value from the group with the higher average FPKM is higher 
#than the highest value in the other group. In other words, it determines if the two groups/conditions have overlapping 
#ranges. If the ranges overlap, samples are removed from either group until an optimal IGSI is found or a cut-off threshold 
#is breached and no IGSI is calculated. 

#Part 3. IQR Outlier Analysis
#Outliers are identified using Interquartile Range (IQR) (Tukey 1977). Samples are examined on a per-gene basis to
#identify moderate outliers (vale = IQR x 1.5-3.0) and extreme outliers (value > IQR x 3.0) (Barbato 2010). To reduce skewness
#of the data, the samples are log2 transformed. Prior to transformation, a value of 1 is added to all samples to avoid taking a log of zero. 

#Part 4. Rank-Order Analysis
#The overall rank-order of samples per gene is calculated and averaged across FPKM values. This helps to identify samples that 
#consistently show higher or lower expression across genes. Additionally, the tallies of samples per rank are calculated, to show 
#the overall spread of the data . The program automatically adjusts for +/- log2(fold change) value reported by Cuffdiff, so the 
#order of the samples is maintained between the groups. When fold-change is positive, the order is lowest to highest. When the 
#fold-change is negative, the order of samples is reversed highest to lowest. 


##References

#Barbato G, Barini EM, Genta G, Levi R. (2011). Features and performance of some outlier detection methods. Journal of Applied Statistics, 38, 2133–49.

#Trapnell, C., Hendrickson, D. G., Sauvageau, M., Goff, L., Rinn, J. L., & Pachter, L. (2013). Differential analysis of gene regulation at transcript resolution with RNA-seq. Nature biotechnology, 31(1), 46-53.

#Tukey, J. W. (1977). Exploratory data analysis. Reading, MA: Addison-Wesley Pub. Co.


################################################################################################################################
#############################################Script for Running of Program######################################################

################################################################################################################################
################################################Required Python Modules#########################################################
import math
import numpy
import sys
import collections
import os.path
import copy

################################################################################################################################
############################################Help Menu###########################################################################

#Help Menu for code
def help_menu_prompt ():
    print ("Help Menu\n")
    print ("To Run Program Please Supply the Following Required Parameters\n")
    print ("Program assumes all files are located in the same directory as program")
    print ("")
    print ("Required Parameters")
    print ("--Group_1_Name")
    print ("--Group_2_Name")
    print ("")
    print ("")
    print ("Optional Parameters")
    print ("--Group_1_Max_Exclusion (default setting = 1)")
    print ("--Group_2_Max_Exclusion (default setting = 1)")
    print ("--Fold_Change_Cutoff (default setting = 1.3)")
    print ("--Statistically_Significant_Data_Only (default setting = true)")
    print ("--Input_File_Location (default setting = local directory)")
    print ("--Output_File_Location (default setting = local directory) ")
    print ("")
    print ("")
    print ("The order of the parameters does not matter, but capitalization, underscores and spelling does matter")
    print ("VERY IMPORTANT when putting pathways of the file \\ at the end is extremely important!!!")
    print ("")
    print ("Definition of Parameters")
    print ("Group_1_Name - Name of the first group of samples analyzed by Cuffdiff aka the ""condition"" (see read_groups.info file)")
    print ("Group_2_Name - Naming of the second group of samples analyzed by Cuffdiff aka the ""condition"" (see read_groups.info file)")
    print ("Group_1_Max_Exclusion - The maximum number of samples to exclude from group 1")
    print ("    cannot exceed number of samples in group")
    print ("Group_2_Max_Exclusion - The maximum number of samples to exclude from group 2")
    print ("    cannot exceed number of samples in group")
    print ("Fold_Change_Cutoff - Threshold for fold change between the two groups")
    print ("Statistically_Significant_Data_Only - filter the Cuffdiff dataset for only significant genes or not")
    print ("    acceptable inputs are 'true' or 'false'")
    print ("Input_File_Location - Location of Cuffdiff data ---REMEMBER forward slash /")
    print ("Output_File_Location - Location to write data ---REMEMBER forward slash /")

    print ("")
    print ("EXAMPLE SUBMISSION IN UNIX")
    print ("")
    print ("python name_of_program \\")
    print ("--Group_1_Name Blue_Chickens \\")
    print ("--Group_2_Name Green_Chickens \\")
    print ("--Group_1_Max_Exclusion 1 \\")
    print ("--Group_2_Max_Exclusion 1 \\")
    print ("--Fold_Change_Cutoff 1.3 \\")
    print ("--Statistically_Significant_Data_Only true \\")
    print ("--Input_File_Location your/favorite/directory/ \\")
    print ("--Output_File_Location your/favorite/output/directory/ \\")
    print ("")
    sys.exit()




####################################################Code to Parse the Parameter File################################################
#Definition splits up a line based on tab
#to return the input value from that line
def split_variable (line):
    key,value = line.split("\t")
    return (value)


#Testing the input from the user to see if its a number (decimals acceptable)
def test_Number_Input(value):
    #converts any floats to numbers 
    #and tests using the python built-in types "isdgit"
    #Note: isdigit does NOT accept decimal places, so input needs to be slightly modified)
    test=value.replace('.','',1).isdigit()
    #examines the results (True/False)
    if test==True:
        result='Pass'
        #print (result)
    else:
        result='Fail'
        #print (result)
    #return the results for the function
    return (result)



#Function retrieves parameter information from a supplied paramter file 
def parsing_input_parameter_file(program_parameters):
    
    #default parameters (will be overridden by user---input)
    Group_1_Max_Exclusion=1
    Group_2_Max_Exclusion=1
    Fold_Change_Cutoff=1.3
    Maximum_Percentile_to_Exclude=0.25
    data_Statistical_Filter='true'
    #Currenting working directory default parameters (same directory as program)
    working_directory= os.getcwd()
    Input_File_Location=working_directory+'\\'
    #global Output_File_Location
    Output_File_Location=working_directory+'/' #for UNIX environment this symbol is required and it works fine in PC submission

    file = open (program_parameters, 'r')

    for line in file:
        if line.startswith("--"):
            line=line.rstrip('\n')
            parsed_parameters=line.split("--")
            for x in range(1, len(parsed_parameters)):
                inputs = parsed_parameters[x].split(" ")
                    
                if inputs[0] == "Group_1_Name":
                    Group_1_Name=inputs[1]     

                elif inputs[0] == "Group_2_Name":
                    Group_2_Name=inputs[1]

                elif inputs[0] == "Group_1_Max_Exclusion":
                    Group_1_Max_Exclusion=inputs[1]

                    #Testing user input (verify numeric value)
                    result = test_Number_Input(Group_1_Max_Exclusion)
                    #If the result passes do this or if the result fails do something else
                    if result =='Pass':
                        pass
                    else:
                        print ("")
                        print ("ERROR Alert!")
                        print ("Please check Group_1_Max_Exclusion value- incorrect input")
                        print ("Incorrect input was: ", Group_1_Max_Exclusion) 
                        sys.exit()

                elif inputs[0] == "Group_2_Max_Exclusion":
                    Group_2_Max_Exclusion=inputs[1]
                    
                    #Testing user input (verify numeric value)
                    result = test_Number_Input(Group_2_Max_Exclusion)
                    #If the result passed do this or if the result fails do that
                    if result =='Pass':
                        pass
                    else:
                        print ("")
                        print ("ERROR Alert!")
                        print ("Please check Group_2_Max_Exclusion value- incorrect input")
                        print ("Incorrect input was: ", Group_2_Max_Exclusion)
                        sys.exit()
                    

                elif inputs[0] == "Fold_Change_Cutoff":
                    Fold_Change_Cutoff=inputs[1]

                    #Testing user input (verify numeric value)
                    result = test_Number_Input(Fold_Change_Cutoff)
                    #If the result passes do this or if the result fails do that
                    if result =='Pass':
                        pass
                    else:
                        print ("")
                        print ("ERROR Alert!")
                        print ("Please check Fold_Change_Cutoff value- incorrect input")
                        print ("Incorrect input was: ", Fold_Change_Cutoff)
                        sys.exit()

                elif inputs[0] == "Input_File_Location":
                    Input_File_Location=inputs[1]

                elif inputs[0] == "Output_File_Location":
                    #print("hello")
                    #print(Output_File_Location)
                    Output_File_Location=inputs[1]

                elif inputs[0] == "Statistically_Significant_Data_Only":
                    data_Statistical_Filter=inputs[1]

                elif inputs[0] == "help" or inputs[0]=="h" or inputs[0]=="Help":
                    printing_help_menu=help_menu_prompt()

                ###This ELSE does not really work because of how the parameter file was parsed, mainly
                #built for purposes of using with qsub file or command line input
                else:
                    print (inputs)
                    print ("")
                    print ("ERROR Alert!")
                    print ("Please double check your input parameters, something was not quite right")
                    print ("Type: --help to see a list of options and acceptable input for the program")
                    sys.exit()
        else:
            pass

    return{'Group_1_Name':Group_1_Name, \
           'Group_2_Name':Group_2_Name, 'Group_1_Max_Exclusion':Group_1_Max_Exclusion, \
           'Group_2_Max_Exclusion':Group_2_Max_Exclusion, 'Fold_Change_Cutoff':Fold_Change_Cutoff, \
           'input_File_Location':Input_File_Location, 'Output_File_Location':Output_File_Location, \
           'data_Statistical_Filter':data_Statistical_Filter}


#######################Function to parse apart user input inputed on the qsub file or command line##############################
#function parses apart a user input from UNIX environment or command line input
def parsing_input(program_parameters):

    #default parameters (will be overridden by user---input)
    Group_1_Max_Exclusion=1
    Group_2_Max_Exclusion=1
    Fold_Change_Cutoff=1.3
    data_Statistical_Filter='true'
    #Currenting working directory default parameters (same directory as program)
    working_directory= os.getcwd()
    Input_File_Location=working_directory+'\\'
    #global Output_File_Location
    Output_File_Location=working_directory+'/' #for UNIX environment this symbol is required and it works fine in PC

    parsed_parameters=program_parameters.split("--")
    #print (parsed_parameters)
    for x in range(1, len(parsed_parameters)):
        inputs = parsed_parameters[x].split(" ")
        
        if inputs[0] == "Group_1_Name":
            Group_1_Name=inputs[1]     

        elif inputs[0] == "Group_2_Name":
            Group_2_Name=inputs[1]

        elif inputs[0] == "Group_1_Max_Exclusion":
            Group_1_Max_Exclusion=inputs[1]

            #Testing user input (verify numeric value)
            result = test_Number_Input(Group_1_Max_Exclusion)
            #If the result passed do this or if the result fails do that
            if result =='Pass':
                pass
            else:
                print ("")
                print ("ERROR Alert!")
                print ("Please check Group_1_Max_Exclusion value- incorrect input")
                print ("Incorrect input was: ", Group_1_Max_Exclusion)
                sys.exit()

        elif inputs[0] == "Group_2_Max_Exclusion":
            Group_2_Max_Exclusion=inputs[1]

            #Testing user input (verify numeric value)
            result = test_Number_Input(Group_2_Max_Exclusion)
            #If the result passes do this or if the result fails do that
            if result =='Pass':
                pass
            else:
                print ("")
                print ("ERROR Alert!")
                print ("Please check Group_2_Max_Exclusion value- incorrect input")
                print ("Incorrect input was: ", Group_2_Max_Exclusion)
                sys.exit()
                    
        elif inputs[0] == "Fold_Change_Cutoff":
            Fold_Change_Cutoff=inputs[1]

            #Testing user input (verify numeric value)
            result = test_Number_Input(Fold_Change_Cutoff)
            #If the result passes do this or if the result fails do that
            if result =='Pass':
                pass
            else:
                print ("")
                print ("ERROR Alert!")
                print ("Please check Fold_Change_Cutoff value- incorrect input")
                print ("Incorrect input was: ", Fold_Change_Cutoff)
                sys.exit()
                
        elif inputs[0] == "Statistically_Significant_Data_Only":
            data_Statistical_Filter=inputs[1]

        elif inputs[0] == "Input_File_Location":
            Input_File_Location=inputs[1]

        elif inputs[0] == "Output_File_Location":
            Output_File_Location=inputs[1] 

        elif inputs[0] == "help" or inputs[0]=="h" or inputs[0]=="Help":
            printing_help_menu=help_menu_prompt()

        else:
            print ("")
            print ("ERROR Alert!")
            print ("Please double check your input parameters, something was not quite right")
            print ("Type: --help to see a list of options and acceptable input for the program")
            sys.exit()

    return{'Group_1_Name':Group_1_Name, \
           'Group_2_Name':Group_2_Name, 'Group_1_Max_Exclusion':Group_1_Max_Exclusion, \
           'Group_2_Max_Exclusion':Group_2_Max_Exclusion, 'Fold_Change_Cutoff':Fold_Change_Cutoff, \
           'input_File_Location':Input_File_Location, 'Output_File_Location':Output_File_Location, \
           'data_Statistical_Filter':data_Statistical_Filter}

################################################################################################################################
###############################################Part 1. Parsing FPKM Data #######################################################

## Function reads in a file and splits on the tab
def read_input_split(inputfile):
    inputf = [inp.rstrip() for inp in open(inputfile)]
    for i in range(len(inputf)):
        inputf[i] = inputf[i].split('\t')
    return inputf


#Function parses the sample IDs and CuffDiff IDs into variables and lists that are utilized
#by a later function when matching FPKM values to the actual sample ID
#Note: CuffDiff re-labels samples during analysis, so a key needs to be created to re-convert back to an understable ID
def parsing_read_Groups(groupList):
    grp = [] #names of affected and unaffected groups
    idmatch = []
    ididx = []
    for i in range(1,len(groupList)): #length of the glist is the number of lines (skips header)
        if not grp.count(groupList[i][1]): #retrieves variable from line "i" from index 1 (Unaffected versus Affected categories)
            grp.append(groupList[i][1]) 
        cuffid  = groupList[i][1]+groupList[i][2] #Creates Cuffdiff IDs (condition + replicate_num)
        #print (cuffid)
        idmatch.append([]) #creates an empty list in a list
        ididx.append(cuffid) #List of all samples (Cuffdiff ids)
        idmatch[i-1].append(cuffid) #adds to first empty list
        ids = groupList[i][0].split('/')[0] #Splitting out Sample IDs from data
        idmatch[i-1].append(groupList[i][1]+'_'+ids) #Creates a list of lists ("Cuff Diff Ids" and "Samples IDs" example: ['HM_U0', 'HM_U_Sample_435']
    return{'grp':grp, 'idmatch':idmatch, 'ididx':ididx}



#Function takes in CuffDiff data and converts the gene lists and output into a file format that is more
#user friendly and understandable for later parsing.
#Function calls the "parsing_read_Groups" function to create a matching key for the data and then begins parsing
#genes.read_group_tracking file in a series for loops to create an output where list consists of genes (y-axis)
#and corresponding sample FPKM values (x-axis)
#Three input parameters: genes.read_group_tracking (Cuffdiff output), read_groups.info (Cuffdiff output), Output_File_Location)
def cfdirdcksing2(infile,read_Group_File,outfile,Output_File_Location):
    outf = open(Output_File_Location+outfile,'w') #setting up new file with parsed FPKM values
    outf.write('Gene') #adding first title to output file

    #Parsing of the Cuffdiff read_Group_File (IDs)
    grouplist = read_input_split(read_Group_File) #all the data from the read_groups.info read into this variable, can be read line by line
    read_Groups=parsing_read_Groups(grouplist)
    grp=(read_Groups['grp'])
    idmatch=(read_Groups['idmatch'])
    ididx=(read_Groups['ididx'])

    #Begin Parsing of Cuffdiff genes.read.group_tracking file
    inf1 = read_input_split(infile) #Cuffdiff genes.read.group_tracking file
    nogrp =len(grp) #number of groups
    noidx = len(idmatch) #number of samples
    gn1 = [] #list of genes in analysis
    nog1 = 0
    for i in range(1,len(inf1)): #genes_read_group_tracking length of file
        if not gn1.count(inf1[i][0]):
            #First line of the file add in Sample IDs from list of list to file
            if nog1 ==0:
                for x in range(noidx):
                    outf.write('\t'+idmatch[x][1])
                outf.write('\n')
            #prints the samples FPKM to the file in correct order of samples
            else:
                for x in range(noidx):
                    outf.write('\t'+gexp[x]) #filled in gexp list from all samples
                outf.write('\n')
            gexp = [] #setups an empty list to place FPKM values into
            for z in range(noidx):
                gexp.append('NA') #fill in list with NAs, which will be replaced (place markers in case empty values)
            cuffid = inf1[i][1]+inf1[i][2] #retrieves sample ID
            iidx = ididx.index(cuffid)
            #print (i,cuffid,iidx)
            gexp[iidx] = inf1[i][6] #place first value into gexp list based index assignment
            outf.write(inf1[i][0]) #adds the gene name to the file
            gn1.append(inf1[i][0]) #add to gene list
            nog1+=1 #gene list counter
        #retrieves the rest of the samples data for a specific gene
        #from the file (first sample already placed in gexp list in correct index)
        else:
            cuffid = inf1[i][1]+inf1[i][2] #identify sample INDEX
            iidx = ididx.index(cuffid)
            gexp[iidx] = inf1[i][6] #place sample FPKM in correct INDEX location

    #prints the last line of FPKM values to data file---needed to get the last bit of data
    for i in range(noidx):
        outf.write('\t'+gexp[i])
    outf.write('\n')
    outf.close()

################################################################################################################################
#######################################Part 2. Functions to Run Inter-Group Seperation Index (IGSI) Analysis####################

#Function pulls in a list of items and retrieves the 2nd item in the list
def getKey(item):
    return item[1]


## Finds the average of a list of numbers
def average(x):
    assert len(x) > 0
    return float(sum(x)) / len(x)


#Function parses apart the expression data file prior created in the program and 
# returns lists of the samples and their expression values for one gene
def parsing_expression_data (inputf, i, nog1, nog2, g1exp, g2exp, g1val, g2val, g1idx, g2idx):
    for j in range(nog1):
        #creating a list inside a list
        g1exp.append([])
        #adding sample names to the group 1 list
        g1exp[j].append(inputf[0][g1idx[j]])
        #adding the sample expression to the group 1 list [[sample 1, value][][]...]
        g1exp[j].append(float(inputf[i][g1idx[j]]))
        #adding expression values to a seperate group 1 list (values only)
        g1val.append(float(inputf[i][g1idx[j]]))
    for j in range(nog2):
        #creating a list inside a list
        g2exp.append([])
        #adding sample names to the group 2 list
        g2exp[j].append(inputf[0][g2idx[j]])
        #adding the sample expression to the group 2 list [[sample 6, value][][]...]
        g2exp[j].append(float(inputf[i][g2idx[j]]))
        #adding expression values to a seperate group 2 list (values only)
        g2val.append(float(inputf[i][g2idx[j]]))
    return {'group_1_expression':g1exp, 'group_1_values':g1val, 'group_2_expression':g2exp, 'group_2_values':g2val}



# Function uses "user" supplied names of group1 and group2 samples (Cuffdiff conditions) and sorts through the header
# of the file for the samples names that match the the groups and returns the list of the indexes
# returns lists of the two group indexes (g1idx and g2idx), the function also returns the ID list
# from the header of the file
def retrieve_sample_Indexes(g1head, g2head, g1idx, g2idx, idlist, inputf):
    #cycles through the header
    for i in range(len(inputf[0])):
        #retrieves samples that start with group_1 sample names
        if inputf[0][i].startswith(g1head):
            g1idx.append(i)
            idlist.append(inputf[0][i])
        #retrieves samples that start with group_2 sample names
        elif  inputf[0][i].startswith(g2head):
            g2idx.append(i)
            idlist.append(inputf[0][i])
    if len(g1idx) == 0 or len(g2idx) == 0:
        print ("")
        print ("ERROR Alert!")
        print ("Please double check your group 1 and group 2 sample names")
        print ("Program did not identify one of the group sets in your dataset")
        print ("User Supplied Group 1 Name: ", g1head)
        print ("User Supplied Group 2 Name: ", g2head)
        sys.exit()
    else:
        pass
                      
    return{'group_1_index':g1idx, 'group_2_index':g2idx, 'id_list':idlist}


# Input file created from the prior function used to parse Cuffdiff Output files
# Function sorts through the prior output data and prioritizes the gene list from Cuffdiff
# g1head = group 1 header
# g2header = group 2 header
# function loops through to the find the best combination of samples that have largest change in expression
# based on various combinations of removing outliers from the data (top and bottom of both datasets) 
def checkindexp(indexpdata,g1head,g2head,g1exmax,g2exmax,tfc):
    
    inputf = read_input_split(indexpdata)

    #Creating files of output data 
    out_log_file=open(indexpdata[:-4]+'_IGSI_gene_log.txt','w')
    out_file_gene_results= open(indexpdata[:-4]+'_IGSI_FoldChange.txt','w')
    
    #Printing headers to the gene expression list file
    out_file_gene_results.write('Gene\tIGSI_Value\tCoeffVar_Grp_1\tCoeffVar_Grp_2\tFlag_Warning\tIGSI_Fold-Change(Grp1/Grp2)\tIGSI_log2FC(Grp1/Grp2)\tIGSI_Abs(log2FC(Grp1/Grp2))\tGrp1_Exclude\tGrp2_Exclude\n')
    
    #Starting setting up the program to Run
    #Group 1 Index 
    g1idx = []
    #Group 2 Index
    g2idx = []
    #List of Samples IDs
    idlist= []

    #length of the entire file (including header)
    norow = len(inputf)

    #Parses apart header from file 
    sample_indexes=retrieve_sample_Indexes(g1head, g2head, g1idx, g2idx, idlist, inputf)
    #Calls data from dictionary of sample_indexes function results
    g1idx =sample_indexes['group_1_index']
    g2idx =sample_indexes['group_2_index']
    idlist =sample_indexes['id_list']
    
    #count of number of samples in group 1
    nog1 =len(g1idx)
    #count of number of samples in group 2
    nog2 = len(g2idx)

    #Count the number of genes examined
    gene_counter=0
    
    #Print the sample header for log file (all combos possible and also retrieve number of combos)
    numb_of_combos=printing_log_file_header(out_log_file, nog2, nog1, g2head, g1head, g1exmax, g2exmax)

    #Calculate Lowest IGSI value possible (Quick Notes on Lowest IGSI)
    #nog1 is the number of samples of group 1
    #nog2 is the number of samples of group 2
    #g1exmax is the maximum number of samples to exclude group 1(one was added to correct range issues)
    #g2exmax is the maximum number of samples to exclude group 2(one was added to correct range issues)      
    lowest_IGSI=(nog2-(g2exmax-1))/nog2 + (nog1-(g1exmax-1))/nog1

    #Gene No combinations found counter
    gene_no_combo_counter=0
    
    #exclusion individual list
    exind = []
    #retrieves data from the prior created file from gene expression
    for i in range(1,norow):

        #Counters in Dataset
        gene_counter+=1

        #Counts how many times "no combos" is found in analysis
        no_combo_found_counter=0

        #retrieves the gene ID from the file
        out_log_file.write(inputf[i][0])
        
        #group 1 gene expressions list
        g1exp = []
        #group 2 gene expressions list
        g2exp = []
        #group 1 gene expression values ONLY list
        g1val = []
        #group 2 gene expression values ONLY list
        g2val = []
        
        expression_data=parsing_expression_data (inputf, i, nog1, nog2, g1exp, g2exp, g1val, g2val, g1idx, g2idx)
        g1exp=expression_data['group_1_expression']
        g1val=expression_data['group_1_values']
        g2exp=expression_data['group_2_expression']
        g2val=expression_data['group_2_values']

        #Note: Sorted is a built in python function
        #Sorts the group 1 gene list by expression values
        g1exp = sorted(g1exp,key=getKey)
        #Sorts the group 1 gene list by expression values
        g2exp = sorted(g2exp,key=getKey)
        #Sorts the group 1 value list by values 
        g1val = sorted(g1val)
        #print ("The g1val are:", g1val)
        #Sorts the group 1 value list by values
        g2val = sorted(g2val)
        #print ("The g2val are:", g2val)
                
        #group One Mean
        g1mean = average(g1val) 
        #group Two Mean 
        g2mean = average(g2val) 
        
        #Sets initial starting values
        #Max match ratio value
        #Note: Max match ratio (maxmr) is later renamed IGSI value 
        maxmr = 0.0
        #Best fold change
        fold_change = 0.0
        #Group 1 Expression levels
        g1ex = 0
        #Group 2 Expression levels
        g2ex = 0
        #Group 1 Coefficient of Variation
        cv_grp1_var=0
        #Group 2 Coefficient of Variation
        cv_grp2_var=0

        #final flag that gets passed
        final_flag='stuff'
        
        if g1mean >= g2mean:
            #judging values
            judval = 1
        else:
            #judging values
            judval = 2

        #setting up the "for" loop of including/excluding various samples from the top and
        #bottom of the two groups to see if there is a significant difference in values
        for j in range(g2exmax): ###Group 2 Exclusion Maximum
            if judval == 1:
                #counts DOWN expression list, final value is lowest value in list
                if j != 0: #excludes max value
                    g2com = g2val[:-j]
                #captures full group 2 expression list
                else:
                    g2com = g2val
            #group 2 expression higher, so counts up  list, final value is highest value in list
            else: 
                g2com = g2val[j:]
                #print ("The j value is:", j)
                #print ("The g2com values:", g2com)

            for m in range(g1exmax): 
                if judval == 1:
                    #group 2 expression higher, so counts up  list, final value is highest value in list)
                    g1com = g1val[m:]
                else:
                    #counts DOWN expression list, final value is lowest value in list
                    if m != 0: #excludes max value
                        g1com = g1val[:-m]
                        
                    else:
                        #captures full group 1 expression list
                        g1com = g1val 

                #Threshold and Match Ratio (IGSI Testing)
                # Tests various combination of the two groups and if they pass 
                # the threshold where maximum value of one of the groups is less than the minimum value of the other group       
                if (judval == 1 and min(g1com) > max(g2com)) or (judval == 2 and min(g2com) > max(g1com)):
                    #Find the max ratio (IGSI value) for the combination of the two groups
                    #The maximum match ratio possible is 1 + 1 = 2
                    match_test = (nog2-j)/nog2+(nog1-m)/nog1 #Example calculation ((6-3)/6) + ((8-2)/8) = 1.25
                    out_log_file.write('\t'+str(match_test))

                    ###############################################
                    #Fold Change Results (First analyzing if no zeros in average)
                    if average(g2com) != 0.0 and average(g1com) !=0.0:
                        foldchange = average(g1com)/average(g2com)
                        group_1_std=numpy.std(g1com)
                        group_2_std=numpy.std(g2com)

                        #Calculating Co-efficent of variant in Data
                        cv_group_1=group_1_std /(average(g1com))*100
                        cv_group_2=group_2_std /(average(g2com))*100

                        out_log_file.write('|'+str(foldchange))
                        out_log_file.write('|Pass')
                    
                        flag = 'Pass'
                        
                    #Dealing with Zeros in the datasets
                    else:
                        #If group 2 has no expression level, foldchange is recorded with a flag in dataset
                        # Example output: 1.375|FLAG-value|fold_change|Pass
                        if average (g2com) == 0.0:

                            foldchange = average(g1com)/0.001

                            group_1_std=numpy.std(g1com)
                            group_2_std= 0

                            cv_group_1=group_1_std /(average(g1com))*100
                            cv_group_2= 0

                            out_log_file.write('|'+str(foldchange))
                            out_log_file.write('|FLAG-Mean_Grp_2_Zero')

                            flag = 'FLAG-Mean_Grp_2_Zero'
                    
                        #If group 1 has no expression level, foldchange is 0 and its recorded as LFC 
                        #Example output 1.708|0.0|LFC
                        else:
                            foldchange = 0.001/average(g2com)

                            group_1_std=0
                            group_2_std=numpy.std(g2com)

                            cv_group_1=0
                            cv_group_2=group_2_std /(average(g2com))*100

                            out_log_file.write('|'+str(foldchange))
                            out_log_file.write('|FLAG-Mean_Grp_1_Zero')

                            flag = 'FLAG-Mean_Grp_1_Zero'
                    
                    #Comparing the fold change to the cutoff provided by the user input
                    if abs(math.log(foldchange,10))> math.log(tfc,10):
                        out_log_file.write('|Fold_Change_Sig')
                        #finding the group combination with largest match ratio (IGSI) value (greatest value possible =2)
                        if maxmr < match_test:
                            #added another filter level to find best value
                            #that finds the match ratio (IGSI) with greatest foldchange
                            if fold_change < foldchange:
                                maxmr = match_test
                                fold_change = foldchange
                                g1ex = m
                                g2ex = j
                                cv_grp1_var=cv_group_1
                                cv_grp2_var=cv_group_2
                                final_flag=flag
                    else:
                        out_log_file.write('|Fold_Change_Not_Sig') 
                else:
                    out_log_file.write('\tNo_Combination')
                    no_combo_found_counter+=1
                    #print ("No_Match")
        if no_combo_found_counter==numb_of_combos:
            gene_no_combo_counter+=1


        out_log_file.write('\n')
    
        sample_exclusion_list=print_gene_results(inputf, i, out_file_gene_results, maxmr, fold_change, g1ex,
                                                 g2ex, judval, g1exp, exind, g2exp, cv_grp1_var, cv_grp2_var, final_flag, lowest_IGSI)

    out_log_file.close()
    out_file_gene_results.close()
    
    #Printing Exclusion List
    printing_Gene_Exclusion_Per_Sample(indexpdata, idlist, sample_exclusion_list, gene_counter, gene_no_combo_counter)

#function prints the header to log file for the two groups being compared
def printing_log_file_header(out_log_file, nog2, nog1, g2head, g1head, g1exmax, g2exmax):
    numb_of_combos=0
    #prints the header
    out_log_file.write('Gene')
    #starts the range of j based on number in group 2
    for j in range(g2exmax):
        #starts the range of m based on number in group 1
        for m in range(g1exmax):
            #prints out the various combination of genes 
            out_log_file.write('\t'+g1head+str(nog1-m)+'+'+g2head+str(nog2-j))
            #make a list of combinations 
            numb_of_combos+=1
    out_log_file.write('\n')

    return (numb_of_combos)

#function prints the total number of genes excluded per sample, by accepting the list
# of samples and then counting the number of times that sample appears in the list
def printing_Gene_Exclusion_Per_Sample(indexpdata, idlist, sample_exclusion_list, gene_counter, gene_no_combo_counter):
    
    #Creating Output file Number of Genes Excluded Per Sample
    #out_file_ind=open(indexpdata[:-4]+'_SUMMARY_REPORT_FILE.txt','w')
    out_file_ind=open(indexpdata[:-4]+'_summary_report.txt','w')

    out_file_ind.write('SUMMARY REPORT OF ANALYSIS\n')
    out_file_ind.write('\n')

    #Writing headers of Two of the files
    out_file_ind.write('No. of times samples were excluded from optimal grouping\n')
    out_file_ind.write('Sample_ID\tCounts_of_Exclusion\tNo_Exclusion_Count(%)\n')
    
    for i in range(len(idlist)):
        #prints the sample name
        out_file_ind.write(idlist[i])
        #counts the number of times a sample shows up in the list of exclusion
        noindex = sample_exclusion_list.count(idlist[i])
        #prints the excluded sample count
        out_file_ind.write('\t'+str(noindex)+'\t'+'{:.2%}'.format(noindex/(gene_counter-gene_no_combo_counter)))
        out_file_ind.write('\n')

    out_file_ind.write('The number of genes analyzed was: '+str(gene_counter)+'\n')
    out_file_ind.write('The number of genes where no IGSI value could be calculated: '+str(gene_no_combo_counter)+'\n')
    out_file_ind.write('The number of genes with a passing IGSI value: '+str(gene_counter-gene_no_combo_counter))   
    out_file_ind.write('\n')    
    out_file_ind.write('\n')
    out_file_ind.write('\n')
    out_file_ind.close()

#function prints out the gene expression results for comparison of samples (affected versus unaffected)
def print_gene_results(inputf, i, out_file_gene_results, maxmr, fold_change, g1ex, g2ex, judval,
                       g1exp, exind, g2exp, cv_grp1_var, cv_grp2_var, final_flag, lowest_IGSI):
    #prints the gene name to file
    out_file_gene_results.write(inputf[i][0])

    #dealing with No Combination Found in the data situations
    if maxmr == 0.0:
        out_file_gene_results.write('\t'+"<"+str(lowest_IGSI)+'\tNA\tNA\tNo_Combination\t'
                                    'NA\tNA\tNA\tNA\tNA\n')
    else:
        #prints the max matching test value to file  
        out_file_gene_results.write('\t'+str(maxmr))

        #prints the coefficient of variation for group 1 and group 2
        out_file_gene_results.write('\t'+str(cv_grp1_var)+"%")
        out_file_gene_results.write('\t'+str(cv_grp2_var)+"%")

        out_file_gene_results.write('\t'+final_flag)
        
        #prints the best fold change value
        out_file_gene_results.write('\t'+str(fold_change))
        if fold_change != 0.0:
            out_file_gene_results.write('\t'+str(math.log(fold_change,2))+'\t'+str(abs(math.log(fold_change,2))))
        else:
            out_file_gene_results.write('\t.')

        #Prints the range of samples that were excluded from analysis
        #The printing depends on the "judging value" which was determined based on which
        #grouping had the higher average expression level
        #Now the exind list is collecting all the samples excluded for future counting

        #Prints the number of samples excluded for group 1
        out_file_gene_results.write('\t'+str(g1ex))    
        for m in range(g1ex):
            if judval == 1:
                out_file_gene_results.write(','+g1exp[m][0])
                exind.append(g1exp[m][0])
            else:
                out_file_gene_results.write(','+g1exp[-m-1][0])
                exind.append(g1exp[-m-1][0])
        #Prints the number of samples excluded for group 2
        out_file_gene_results.write('\t'+str(g2ex))
        for m in range(g2ex):
            if judval == 1:
                out_file_gene_results.write(','+g2exp[-m-1][0])
                exind.append(g2exp[-m-1][0])
            else:
                out_file_gene_results.write(','+g2exp[m][0])
                exind.append(g2exp[m][0])

        
        
        out_file_gene_results.write('\n')

    return (exind) #return list of excluded samples for sample exclusion output

################################################################################################################################
#########################################################Functions to Merge and Filter Datasets#################################

#Parsing through file for header and data
#Function accepts a file and starts to parse through the data
#returning the results as a retrievable dictionary (header and dictionary of lines
#first value in the line is the key and the entire line is the value
#Dictionary of lines can be used to compare common lines between files
def retrieve_data (file_name):
    file= open(file_name, 'r')
    file_dic={}
    for line in file:
        if line.startswith("Gene"):
            line=line.rstrip('\n')
            header_info=(line)
        elif line.startswith("test_id"):
            line=line.rstrip('\n')
            header_info=(line)
        else:
            line=line.rstrip('\n')
            values = line.split("\t")
            key = (values[0])
            values=line
            file_dic.update({key:values})
    return{'header_info':header_info, 'file_dic':file_dic}


#Parsing through file for header and data---but REMOVES first columns from data
#Function accepts a file and starts to parse through the data
#returning the results as a retrievable dictionary (header and dictionary of lines
#first value in the line is the key and the entire line is the value
#Dictionary of lines can be used to compare common lines between files
def retrieve_data_modified (file_name):
    file= open(file_name, 'r')
    file_dic={}
    for line in file:
        if line.startswith("Gene"):
            line=line.rstrip('\n')
            header = line.split("\t")
            header=header[1:]
            modified_header=('\t'.join(map(str,header)))
    
        else:
            line=line.rstrip('\n')
            values = line.split("\t")
            key = (values[0])
            values=values[1:]
            modified_values=('\t'.join(map(str,values)))   
            file_dic.update({key:modified_values})
    return{'header_info':modified_header, 'file_dic':file_dic}


#Matching the Merged Data Cuff Diffs "gene_exp.diff" file and the output from
#sorted FPKM data. The data gets sorted as a dictionary (key=gene) and then the two files
#are compared for overlapping based on the statistical significance setting. If setting is "true"
# only significant genes from CuffDiff are outputed from the FPKM data for further analysis,
# if the setting is "false" all genes are analyzed. 
def filter_Cuff_Diff(cuffDiff_Gene_Stats_Reports, statistical_filter, statistical_Name,Output_File_Location):
    #parsed_data_file=Output_File_Location+statistical_Name+'.txt'
    out_file_Filtered=open(Output_File_Location+statistical_Name+'.txt','w')

    #Getting data from the parsed FPKM data file
    parsed_gene_Expression_Data=retrieve_data(Output_File_Location+"FPKM_Samples.txt")
    parsed_Gene_Header=parsed_gene_Expression_Data['header_info']
    parsed_Gene_Dic=parsed_gene_Expression_Data['file_dic']

    #Getting data from the gene_exp.diff from CuffDiff
    gene_Exp_Diff=retrieve_data(cuffDiff_Gene_Stats_Reports)
    gene_Exp_Diff_Header=gene_Exp_Diff['header_info']
    gene_Exp_Dic=gene_Exp_Diff['file_dic']

    out_file_Filtered.write(parsed_Gene_Header +'\n')
    #Nested loops to sort through the data, but only prints one files results based on the other files data
    for x in (gene_Exp_Dic):
        #Values are not tab seperated and need to parsed to find the correct column
        gene_Exp_Data = (gene_Exp_Dic[x].split("\t"))
        #Column (index 13) is where a gene is reported as significant, identifying results and filtering
        # all non-significant data values
        if gene_Exp_Data[13]==statistical_filter:
            continue
            print (gene_Exp_Dic[x])
        #identifying genes found in the parsed FPKM file that overlap with significant genes and only
        # printing out significant genes to be analyzed in the program
        for y in (parsed_Gene_Dic):
            if x == y :
                out_file_Filtered.write(parsed_Gene_Dic[y] + '\n')
            else:
                pass 

    out_file_Filtered.close()


#Creating Output file for Merged Data Cuff Diffs "gene_exp.diff" file and the output from
#Match Test (IGSI Test) output. The data gets sorted as a dictionary (key=gene) and then the two files
#are merged together in one giant file. Only matching data is outputed, so "non-significant" genes will get
#filtered out
def merge_Cuff_Diff_Match(cuffDiff_Gene_Stats_Reports, statistical_Name,Output_File_Location):
    
    out_file_Merge=open(Output_File_Location+statistical_Name+'_gene_exp_IGSI.txt','w')

    #Getting data from the Match_Test file
    match_File_Data=retrieve_data_modified(Output_File_Location+statistical_Name+'_IGSI_FoldChange.txt')
    match_File_Header=match_File_Data['header_info']
    match_File_Dic=match_File_Data['file_dic']

    #Getting data from the gene_exp.diff from CuffDiff
    gene_Exp_Diff=retrieve_data(cuffDiff_Gene_Stats_Reports)
    gene_Exp_Diff_Header=gene_Exp_Diff['header_info']
    gene_exp_dic=gene_Exp_Diff['file_dic']

    out_file_Merge.write(gene_Exp_Diff_Header+'\t'+match_File_Header+'\n')
    #Nested loops to sort through and merge the final datasets where the genes overlap
    for x in (gene_exp_dic):
        for y in (match_File_Dic):
            if x == y :
                out_file_Merge.write(gene_exp_dic[x] + '\t' + match_File_Dic[y] + '\n') 
            else:
                pass 

    out_file_Merge.close()

####################################################################################################################################################
#################################################Part 3. Functions to Run Outlier Analysis of the Data##############################################
    
# Outlier function to detect the outlier datapoints
# The function follows the IQR (Interquartile Range Method) developed by Tukey 1977 (Exploratory Data Analysis) 
# we classified the outliers as moderate or extreme based on the criteria from (Barbato 2010)
# Link for paper: http://www.tandfonline.com/doi/abs/10.1080/02664763.2010.545119
# It should be recognized the python code was adapted from the basic
# setup found on Colin Gorrie's blog http://colingorrie.github.io/outlier-detection.html
# which takes advantage of numpy's percentile module, but greatly modified and expanded
# for our specific data input
def outliers_iqr_test(values):
    #Note numpy calculates percentiles slightly different(weighted), need to add 'midpoint' to calculate IQR without weight correction
    quartile_one=numpy.percentile(values, 25, interpolation = 'midpoint')
    
    quartile_third=numpy.percentile(values, 75, interpolation = 'midpoint')
    
    iqr = quartile_third - quartile_one

    #moderate ranges for whiskers
    upper_range = quartile_third + (1.5*iqr)
    lower_range = quartile_one - (1.5*iqr)

    #extreme ranges for whiskers
    extreme_upper_range=quartile_third + (3.0*iqr)
    extreme_lower_range = quartile_one - (3.0*iqr)
    
    #return a dictionary of variables for testing/printing
    return {'quartile_one':quartile_one, 'quartile_third':quartile_third,
            'lower_range':lower_range, 'upper_range':upper_range,
            'extreme_upper_range':extreme_upper_range, 'extreme_lower_range':extreme_lower_range}

# function opens up the data file and retrieves the list of
# sampels from the header of the file and creates an empty
# sample dictionary to be filled in at a later point
def sample_outlier_dic_func(file_name):
    file= open(file_name, 'r')
    sample_outlier_dic={}
    for line in file:
        if line.startswith("Gene"):
            line=line.rstrip('\n')
            line=line.split("\t")
            samples=line[1:]
            for x in range(0, len(samples)):
                sample_outlier_dic.update({samples[x]:[0,0]})
        else:
            pass
    return (sample_outlier_dic)

#function retrieves the genes and their corresponding IGSI value from
# output file and creates a dictionary of this information
def retrieve_igsi_data(file_name):
    #gene:igsi dictionary 
    igsi_dict={}
    file= open(file_name, 'r')
    for line in file:
        if line.startswith("Gene"):
            continue
        else:
            line=line.rstrip('\n')
            values = line.split("\t")
            key = (values[0])
            igsi_value=values[1]
            igsi_dict.update({key:igsi_value})
    return(igsi_dict)


# function pulls in the data file and empty sample dictionary and begins
# analyzing all the samples for outlier data points. When a point is
# found a count gets added to that samples outlier counter found in
# the outlier sample dictionary
def outlier_test(parsed_data_file, sample_outlier_dic, Output_File_Location, statistical_Name, igsi_dict ):
    #Creating a log file to validate data
    out_outlier_Log=open(Output_File_Location+statistical_Name+'_Outlier_log.txt','w')
    out_outlier_Log.write("Gene\tIGSI_value\tLog2(Data_Values+1)\tFirst_Quartile\tThird_Quartile\tWhisker_Range\tExtreme_Whisker_Range\tOutlier_Values\n")
    file= open(parsed_data_file, 'r')

    #count the number of genes analyzed
    gene_counter=0

    for line in file:
        #flag for printing outlier status of data
        outlier_flag=0
        #Begin sorting through file, identify header with sample name
        if line.startswith("Gene"):
            line=line.rstrip('\n')
            line=line.split("\t")
            #retrives a list of the samples from the data
            samples=line[1:]
            continue
        #starting to parse through the FPKM values of the data for outliers 
        else:
            #progress counter
            gene_counter+=1

            #start parsing the gene FPKM data
            line=line.rstrip('\n')
            line=line.split("\t")
            #retrieves the gene name
            gene=line[0]
            #retrives the FPKM values (index same as index of samples---easy to 
            data=line[1:]
            #creating a list of the values for record keeping
            data_values=','.join(data)

            igsi_value=igsi_dict[gene]

            #converts the string of values to list of float numbers
            #However first add "1" to ALL the values to deal with zeroes in the dataset for next step.
            #then log transforms the data values for analysis (log transformation to try to normalize data) 
            for x in range(len(data)):
                data[x]=numpy.log2(float(data[x])+1)
                
            #retrieves a list of the IQR range values
            iqr_test_information=outliers_iqr_test(data)
            quartile_one=iqr_test_information['quartile_one']
            quartile_third=iqr_test_information['quartile_third']
            lower_range=iqr_test_information['lower_range']
            upper_range=iqr_test_information['upper_range']
            extreme_upper_range=iqr_test_information['extreme_upper_range']
            extreme_lower_range=iqr_test_information['extreme_lower_range']


            #Joining the log transformed data for printing
            log_transformed_data=[]
            for x in range(len(data)):
                log_transformed_data.append(str(data[x]))
            log_data_values=','.join(log_transformed_data)


            #writes the gene name, a list of the values being analyzed, lower IQR range, upper IQR range and excluded values
            out_outlier_Log.write(gene+'\t'+igsi_value+'\t'+log_data_values+'\t'+str(quartile_one)
                                  +'\t'+str(quartile_third)+'\t'+"["+str(lower_range)+"-"
                                  +str(upper_range)+"]"+'\t'+"["+str(extreme_lower_range)+"-"
                                  +str(extreme_upper_range)+"]"+'\t')

            #searches the list of values to identify outlier values
            #if a outlier is found, it is written to the log file
            for x in range(len(data)):
                #first find possible extreme outliers
                if data[x]>extreme_upper_range or data[x]<extreme_lower_range:
                    key=samples[x]
                    sample_outlier_dic[key][1]+=1
                    out_outlier_Log.write('Extreme-'+key+':'+str(data[x])+' ')
                    #flag used when printing data to log file
                    outlier_flag+=1
                    pass
                #next find posible moderate outliers
                elif data[x]>upper_range or data[x]<lower_range:
                    key=samples[x]
                    sample_outlier_dic[key][0]+=1
                    out_outlier_Log.write('Mild-'+key+':'+str(data[x])+' ')
                    #flag used when printing data to log file
                    outlier_flag+=1
                else:
                    pass
        if outlier_flag>=1:
            out_outlier_Log.write('\n')
        else:
            out_outlier_Log.write('No_Outlier_Data_Identified\n')
   
    return {'sample_outlier_dic':sample_outlier_dic, 'gene_counter':gene_counter}
    out_outlier_Log.close()

# function prints the counts of the outliers, which are being stored
# as dictionary of key, value pairs
def printing_outlier_dictionary(parsed_data_file, sample_outlier_dic, gene_counter): 
    out_outlier_dic=open(parsed_data_file[:-4]+'_summary_report.txt','a')
    out_outlier_dic.write('Values were classified as Outliers based on IQR Outlier Method\n')
    out_outlier_dic.write('To calculate values for outlier testing: value_tested=log2(FPKM+1)\n')
    out_outlier_dic.write('Mild Outliers are values between (1.5-3.0 x IQR)\n')
    out_outlier_dic.write('Extreme Outliers are values > 3.0 x IQR\n')
    out_outlier_dic.write('Outlier Counter \n')
    out_outlier_dic.write('Sample\tModerate_Outliers\tExtreme_Outliers\tTotal_Counts\tTotal_Counts(%)\n')

    #loops through the dictionary to print all the values
    for key, value in (sample_outlier_dic.items()): 
       out_outlier_dic.write(str(key)+"\t"+str(value[0])+"\t"+str(value[1])+"\t"+str(sum(value))+"\t"+'{:.2%}'.format(sum(value)/gene_counter)+"\n") 

    out_outlier_dic.write("The number of genes analyzed was: "+str(gene_counter))
    out_outlier_dic.write("\n")
    out_outlier_dic.write("\n")
    out_outlier_dic.write("\n")
    
    out_outlier_dic.close()


###################################################################################################################################################
##################################################Part 4. Functions to Run Rank Order Analysis#####################################################


#function retrieves the genes and original log2(fold_change) produced by CuffDiff as a key and value
#dictionary which will be utilzied in another part of the program. The goal is specifically
# to know the direction of the fold change, so the order of ranking can be adjusted later
def retrieve_gene_FC_dic(Output_File_Location, statistical_Name):
    gene_FC_dic={}
    file = open(Output_File_Location+statistical_Name+'_gene_exp_IGSI.txt', 'r')
    for line in file:
        if line.startswith("test_id"):
            continue
        else:
            line=line.split("\t")
            gene_FC_dic.update({line[0]:float(line[9])})
    return(gene_FC_dic)

#creating the a dictionary of samples with an empty list of values that will be added to based on their
#sample ranking, the code reads in the "Parsed_Genes_Only.txt" file because it is smaller and easier to
#manipulate faster then the raw parsed FPKM values
def create_Rank_Order_Dict(Output_File_Location, statistical_Name):
    rank_order_Dict={}
    file = open(Output_File_Location+statistical_Name+'.txt', 'r') 
    #retrieving the list of samples from the file 
    for line in file:
        if line.startswith("Gene"):
            line=line.rstrip()
            line=line.split("\t")
            sample_names=line[1:]
        #now just skip everything else
        else:
            continue

        #Creating rank_order_dictionary (key=sample_name, values=list of ranks)
        for x in range(len(sample_names)):
            rank_order_Dict.update({sample_names[x]:[]})
        return{'rank_order_Dict':rank_order_Dict, 'sample_names':sample_names}

#using OrderedDict
# https://docs.python.org/3/library/collections.html#ordereddict-examples-and-recipes
def tally_rank_order(Output_File_Location, statistical_Name, rank_order_Dict, gene_FC_dic, sample_names):

    out_rank_Log=open(Output_File_Location+statistical_Name+'_rank_order_log.txt','w')
    
    file = open(Output_File_Location+statistical_Name+'.txt', 'r') 

    #Creating a list of ranks for output into the log file
    list_of_ranks=[]
    for x in range(len(sample_names)):
        value=str(x+1)
        list_of_ranks.append("rank_"+value)

    #used a mapping join function to make the joining loop take up less code (python shortcuts)
    string_list_ranks='\t'.join(map(str,list_of_ranks))

    #Printing the header of the log file for ranking orders           
    out_rank_Log.write("Gene\tlog2(Fold_Change)\tRank_Order\t"+string_list_ranks+"\n")

    #retrieving the list of samples from the file 
    for line in file:
        gene_dict={}
        ordered_Gene_Dict={}
        
        if line.startswith("Gene"):
            line=line.rstrip()
            line=line.split("\t")
            sample_names=line[1:]
            continue
        else:
            line=line.rstrip()
            line=line.split("\t")
            gene_name=line[0]
            fpkm_values=line[1:]

        #print gene name to file
        out_rank_Log.write(str(gene_name+"\t"))

        #add values to the sample (key) : FPKM (value) dictionary per line of the file
        for x in range(len(sample_names)):
            gene_dict.update({sample_names[x]:float(fpkm_values[x])})

        #printing the Fold_Change to the log file
        out_rank_Log.write(str(gene_FC_dic[gene_name])+"\t")

        #Determining the change of direction for a gene, to try to keep the influence
        #of affected/unaffected the same if log2(FC) positive sort lowest to highest,
        #if log2(FC) negative reverse order sort the list and find highest to lowest
        if gene_FC_dic[gene_name]>0:
            #Writing to file genes rank ordered verdict
            out_rank_Log.write("Forward_Order\t")
            
            #Sorts the dictionary using Collections Module in Python
            ordered_Gene_Dict=collections.OrderedDict(sorted(gene_dict.items(), key=lambda t: t[1]))

        elif gene_FC_dic[gene_name]<0:
            
            #Writing to file if genes rank ordered got flipped
            out_rank_Log.write("Reverse_Order\t")
            
            #Reverses the order of the dictionary (all needed to add is reverse=True statement)
            ordered_Gene_Dict=collections.OrderedDict(sorted(gene_dict.items(), key=lambda t: t[1], reverse=True))

        else:
            out_rank_Log.write("No_Analaysis_Performed\n")
            continue
        
        #Converts the dictionary to a list (dictionaries cannot be indexed)
        list_Ordered=list(ordered_Gene_Dict.items())
    
        #Sorts through the list of samples and adds that "Rank" to sample rank order Dict
    
        for x in range(0, len(list_Ordered)):
            out_rank_Log.write(str(list_Ordered[x][0])+":"+(str(list_Ordered[x][1]))+"\t")
            rank_order_Dict[list_Ordered[x][0]].append(x+1)
        #Need to move log file to a newline for next gene
        out_rank_Log.write("\n")

    out_rank_Log.close()
    return(rank_order_Dict)

####Functions to Run Creating Tally Per Rank of Each Sample###

#splits the variables from ranking log file (sample_name: value)
def split_variable(variable):
    variable=variable.split(':')
    variable=variable[0]
    return (variable)

#create an empty tallying dictionary for samples
def creating_tallying_dict(sample_names, ranks):
        tallying_dict={}
        ranks_list=[]

        #Create an empty list (length of rankings in study)
        #Note because name of list is the same BECAREFUL
        for x in range(len(ranks)):
            (ranks_list.append(0))

        #Loops through and adds list to samples as a dictionary
        #Added the command copy.deepcopy command because each list would get updated
        #at the same time because Python assumed it was the same list, deepcopy makes each copy independent
        for x in range(len(sample_names)):            
            tallying_dict.update({sample_names[x]:copy.deepcopy(ranks_list)})
            
        return (tallying_dict)
        

#Tallies the number of times a sample appears at a certain position in the ranking order
#Ranking allows for easier production of graphs showing how samples perform globablly
def tallying_samples_per_rank(sample_names, Output_File_Location, statistical_Name):
    out_rank_Log=open(Output_File_Location+statistical_Name+'_rank_order_log.txt','r')
    
    for line in out_rank_Log:
        if line.startswith("Gene"):
            line=line.rstrip()
            line=line.split("\t")
            ranks=line[3:]
            tallying_dict=creating_tallying_dict(sample_names,ranks)          
            #print (tallying_dict)

        else:
            line=line.rstrip()
            line=line.split("\t")
            gene_name=line[0]
            fpkm_values=line[3:]
            
            for x in range(len(fpkm_values)):
                sample_name=split_variable(fpkm_values[x])
                tallying_dict[sample_name][x]+=1

    
    return{'tallying_dict':tallying_dict,'ranks':ranks}   

#function prints out the rank order dictionary of the samples and the tallied rank order dicitonary
#appends to the Summary Report File 
def printing_rank_order(parsed_data_file, rank_order_Dict, tallying_dict, ranks):
    
    out_Final_Rank=open(parsed_data_file[:-4]+'_summary_report.txt','a')
    out_Final_Rank.write('\n')
    out_Final_Rank.write('\n')
    out_Final_Rank.write('Overall Summary Report of Average Rank of Samples')
    out_Final_Rank.write('\n')
    out_Final_Rank.write("Sample_Name\tAverage_Rank_Order\tTotal_Genes_Examined\n")

    #loops through the dictionary to print all the values
    for key, value in (rank_order_Dict.items()):
       out_Final_Rank.write(str(key)+"\t"+str(numpy.mean(value))+"\t"+str(len(value))+"\n")

    ###Printing Tallying Ranks in Summary Report File
    out_Final_Rank.write('\n')
    out_Final_Rank.write('\n')
    out_Final_Rank.write('\n')
    out_Final_Rank.write('Tallied Rank Order Of Each Sample Broken Down by Rank in the Dataset\n')
    out_Final_Rank.write("Sample_Name\t")
    for x in range (len(ranks)):
        out_Final_Rank.write(ranks[x]+"\t")
    out_Final_Rank.write("Total")
    out_Final_Rank.write('\n')
    for key, values in (tallying_dict.items()):
       out_Final_Rank.write(str(key)+"\t")
       sum_of_values=[]
       for x in range (len(values)):
           out_Final_Rank.write(str(values[x])+"\t")
           sum_of_values.append(values[x])
       out_Final_Rank.write(str(sum(sum_of_values)))
       out_Final_Rank.write("\n")

    out_Final_Rank.close
    return()
    
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
  
def main():
    
    print ("")
    print ("###########Analysis of RNA-Seq Data Behavior Tool (ARB)###########")

    #User can now supply parameter file or submit using a qsub/sbatch file on HPC environment

    #Program first test to see if the parameter file is present
    verdict=os.path.isfile('Analysis_CD_Data_Parameter_File.txt')

    #If Program finds the parameter file it then parses it
    if verdict == True:
        #Code to pull in data from parameter file (TESTING ONLY)
        parameter_stuff = parsing_input_parameter_file('Analysis_CD_Data_Parameter_File.txt')
        print ("")
        print ("Loading input from User Provided Parameter File")
        print ("")
        
    ###Else code pulls in data from the qsub file/command line
    else: 
        user_input=sys.argv[1:]
        user_input=' '.join(user_input)
        print ("")
        print ("Loading inpug from User Supplied Parameters from Submission File")
        print ("")
        print (user_input)
        parameter_stuff=parsing_input(user_input)

    #Retrieving parameters from data (either inputs accepted)
    sample_group_1_name = parameter_stuff['Group_1_Name']
    sample_group_2_name = parameter_stuff['Group_2_Name']

    group_1_exclusion_max = int(parameter_stuff['Group_1_Max_Exclusion'])+1
    group_2_exclusion_max = int(parameter_stuff['Group_2_Max_Exclusion'])+1
    fold_change_cutoff = float(parameter_stuff['Fold_Change_Cutoff'])
    global input_File_Location
    input_File_Location = parameter_stuff['input_File_Location']
    global Output_File_Location
    Output_File_Location = parameter_stuff['Output_File_Location']
    data_Statistical_Filter=parameter_stuff['data_Statistical_Filter']

    #Testing Print of all variables
    print ("FOLLOWING PROGRAM PARAMETERS")
    print ("The name of sample group 1 is:", sample_group_1_name)
    print ("The name of sample group 2 is:", sample_group_2_name)
    #Fixing the +1 for implementing prior code to report back proper number from user
    print ("The Maximum Number of Samples to Exclude from Group 1 is:", (group_1_exclusion_max-1))
    print ("The Maximum Number of Samples to Exclude from Group 2 is:", (group_2_exclusion_max-1))
    print ("The fold change cutoff for the data is:", fold_change_cutoff)
    print ("The data will be filtered for statistical significance:", data_Statistical_Filter)
    print ("")
    print ("Directory Locations")
    print ("The location of the data files are:", input_File_Location)
    
    print ("The selected output location for the data was:", Output_File_Location)

    #Determining the statistical filter of the dataset
    if data_Statistical_Filter=='true':
        statistical_filter="no"
        statistical_Name="Sig"
        
    if data_Statistical_Filter=='false':
        statistical_filter=''
        statistical_Name="All"

    #If user input for statistical filter is wrong, kills the program
    while data_Statistical_Filter not in ['true', 'false']:
        print ("")
        print ("ERROR Alert!")
        print ("Please check Statistically_Significant_Data_Only value- incorrect input")
        print ("Incorrect input was: ",data_Statistical_Filter)
        print ("")
        print ("Acceptable inputs are 'true' or 'false'")
        sys.exit()

    print ("")
    print ("ARB is now analyzing the data")

    #Input File Locations
    gene_read_group_tracking=input_File_Location+'genes.read_group_tracking'
    read_groups_info=input_File_Location+'read_groups.info'
    cuffDiff_Gene_Stats_Reports=input_File_Location+'gene_exp.diff'

    #Calling Functions for data analysis
    cfdirdcksing2(gene_read_group_tracking, read_groups_info, "FPKM_Samples.txt",Output_File_Location)

    #Merging Parsed_Data with CuffDiff Gene Stats Report of the data based on statistical test setting
    merge_Final_Results=filter_Cuff_Diff(cuffDiff_Gene_Stats_Reports, statistical_filter, statistical_Name,Output_File_Location)

    parsed_data_file=Output_File_Location+statistical_Name+'.txt'
    checkindexp(parsed_data_file,sample_group_1_name,sample_group_2_name,group_1_exclusion_max, group_2_exclusion_max,fold_change_cutoff)

    #Merging of the final results of the data
    merge_Final_Results=merge_Cuff_Diff_Match(cuffDiff_Gene_Stats_Reports, statistical_Name,Output_File_Location)

    #Running All the Outlier Functions of the Code for Analysis
    sample_outlier_dic=sample_outlier_dic_func(parsed_data_file)
    igsi_dict=retrieve_igsi_data(parsed_data_file[:-4]+'_IGSI_FoldChange.txt')
                                 
    sample_outlier_results=outlier_test(parsed_data_file, sample_outlier_dic, Output_File_Location, statistical_Name, igsi_dict)
    sample_outlier_dic=sample_outlier_results['sample_outlier_dic']
    gene_counter=sample_outlier_results['gene_counter']

    printing_outlier_dictionary(parsed_data_file, sample_outlier_dic, gene_counter)

    #Running All the Rank Order Functions of the Code for Analysis
    gene_FC_dic=retrieve_gene_FC_dic(Output_File_Location, statistical_Name)

    #Getting a dictionary of samples and sample list for future usage
    rank_order_Stuff=create_Rank_Order_Dict(Output_File_Location, statistical_Name)
    rank_order_Dict=rank_order_Stuff['rank_order_Dict']
    sample_names=rank_order_Stuff['sample_names']
 
    rank_order_Dict= tally_rank_order(Output_File_Location, statistical_Name, rank_order_Dict, gene_FC_dic, sample_names)

    ####Getting the Ranks for Each Sample at Each Ranking Position
    tallying_all_ranks=tallying_samples_per_rank(sample_names, Output_File_Location, statistical_Name)
    tallying_dict=tallying_all_ranks['tallying_dict']
    ranks=tallying_all_ranks['ranks']

    #Print results from rank order tests
    printing_rank_order(parsed_data_file, rank_order_Dict, tallying_dict, ranks)

    #Removes the parsed gene file created from Cuffdiff data
    os.remove(Output_File_Location+statistical_Name+'.txt') 

    print ("ARB is done running, have a good day!")

main()

################################################################################################################################



