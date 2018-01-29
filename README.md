# ARB: Analysis of RNA-Seq Data Behavior Tool
***
**NOTE:** Program was submitted to Bioinformatics (ISCB, http://www.iscb.org/) for review.

ARB is a Python-based tool for prioritization of differentially expressed genes, identification of misclassified samples, and ranking of expression severity in differential expression analysis experiments using RNA-seq.  

Code was developed in the Abasht Laboratory at the University of Delaware under the supervision of [Dr. Behnam Abasht](http://canr.udel.edu/faculty/behnam-abasht/).

Please contact Dr. Behnam Abasht (abasht@udel.edu) with any questions or concerns. 


## Why Use ARB?

### 1. Gene Prioritization

Normalized expression data is used to calculate the inter-group separation index (IGSI), ARB's unique gene prioritization metric. IGSI is a measure of the degree of overlap in expression values between to conditions. To calculate IGSI, the program first determines if the two groups/conditions have overlapping ranges. If the ranges overlap, samples are removed from either group until either (1) a subset of samples can be assembled in which expression values between the two groups do not overlap or (2) the exclusion ceiling is breached (exclusion ceiling is set by the user). IGSI theoretically ranges from 0 (all samples removed) to 2 (no samples removed), but it is recommended that the program settings be maintained such that the minimum IGSI is at least 1 (50% of samples removed). IGSI can be used alongside log(fold-change) to prioritize genes more effectively than if fold-change were used alone.

### 2. Identification of Misclassified Samples
Individual samples that are found to cluster with the opposing condition across a large proportion of genes might be misclassified. This is measured by the exclusion rate (i.e. for how many genes does a specific sample need to be excluded to prevent overlapping values between the conditions), which is calculated during IGSI analysis. We recommend further investigation of samples that exhibit an exclusion rate of greather than or equal to 50% because they may be misclassified.


### 3. Rank-Order Analysis
The rank-order analysis assesses each sample's expression level relative to all the other samples for each gene after correcting for directionality of the ranking depending on whether or not a gene is upregulated or downregulated. For studies of disease, this might be used to study disease severity or stage.    


## Installation and Use

Download the program file and place it in your local working directory on a personal computer or server.

There are two possible file input types accepted by the program. The program was originally developed using output from Cuffdiff (Trapnell 2010), but was later expanded to be compatible with the output from any major differential analysis software with the "Custom" `Input_Type` setting.  

### Required input files ("Cuffdiff" setting):

1. read\_groups.info (required to match sample IDs with IDs used by Cuffdiff program) 
* genes.read\_group\_tracking (supplies the FPKM values)
* gene\_exp.diff (required to filter genes reported as statistically significant)

### Required input files ("Custom" setting):

1. Gene\_Exp.diff (differential expressional analysis results in the form of the template file provided on our website) 
* norm\_gene\_counts.txt (normalized gene read counts (FPKM, counts, RPKM, etc.) data in the form of the template file provided on our websites)

Note: Template files can be found with the ARB program files on our website. Open up the template files, located in the template folder and follow the instructions to copy and paste your data into the correct format. 


### System requirements

**python3.6** + the following modules 

* math
* sys
* collections
* os.path
* numpy (note NumPy may need to be installed and does not come standard with python)


### To run the program locally using Python3.6 with a parameter file

Supply the following parameters in a file entitled Analysis\_ARB\_Parameter\_File.txt in the same directory as the python program (just copy and paste the three lines below):


<pre>
Parameter File (space separated)

--Group_1_Name &lt;Name_of_Group_1&gt;

--Group_2_Name &lt;Name_of_Group_2&gt;
</pre>


### To run the program from a HPC environment using a submission file

<pre>
python ARB2.0.0.py \

--Group_1_Name &lt;Name_of_Group_1&gt; \

--Group_2_Name &lt;Name_of_Group_2&gt;
</pre>


### To run the program on the command line

```
python ARB2.0.0.py --Group_1_Name <Name_of_Group_1> --Group_2_Name <Name_of_Group_2> 
```


*NOTE:* When running this program on the command line locally or on a server, beware of issues that may occur if the incorrect version of python is called or if the NumPy module is not installed properly. Some computers may automatically call a prior version of python if the command "python" is given, so a user needs to specifically call python3.6 using full pathways and also needs to install required modules (aka NumPy) to that version. NumPy installed on prior versions of python will NOT work!



Below are the default parameters for the program. To change a default parameter, submit using the following format:

```
--<Parameter> <Variable>
```



Dashed lines and the single space are important when submitting parameters. The order of the parameters does not matter, but capitalization, underscores and spelling do.

*VERY IMPORTANT:* when putting pathways of the files remember to put "/" at the end of the pathway!



### Parameters

All parameters are optional unless otherwise indicated.

`--Input_Type` Changes the type of input files a user can supply (default = Cuffdiff), user changes to "Custom" and supplies required files (see template files for required data and format)

`--Group_1_Name` **REQUIRED** Name of the first group of samples a.k.a. the "first condition" (see read_groups.info file if using Cuffdiff)
 
`--Group_2_Name` **REQUIRED** Name of the second group of samples a.k.a. the "first condition" (see read_groups.info file if using Cuffdiff)

`--Group_1_Max_Exclusion` The maximum number of samples to exclude from group 1, cannot exceed number of samples in group (default = 1)

`--Group_2_Max_Exclusion` The maximum number of samples to exclude from group 2, cannot exceed number of samples in group (default = 1)

`--Statistically_Significant_Data_Only` An option  to filter the dataset for only significant genes where q-value is above 0.05 (default = true)

`--Input_File_Location` Path for location of input files ---REMEMBER forward slash / (default = working directory)

`--Output_File_Location` Path for location to write output data ---REMEMBER forward slash / (default = working directory)




### Example submission using all parameters 
*NOTE:* Assuming submission in HPC environment (aka UNIX), hence the backslashes.


<pre>
python Analysis_Cuff_Diff_Data.Version6.2.0.py \

--Input_Type Custom \

--Group_1_Name Blue_Chickens \

--Group_2_Name Green_Chickens \

--Group_1_Max_Exclusion 1 \

--Group_2_Max_Exclusion 1 \

--Statistically_Significant_Data_Only true \

--Input_File_Location your/favorite/directory/ \

--Output_File_Location your/favorite/output/directory/ \
</pre>



## ARB Output Files 

1. **Sig\_gene\_exp\_IGSI.txt**   
IGSI analysis results merged with original differential expression analysis results for easy sorting. If no IGSI value was found for a gene, it is recorded as "<" the lowest possible IGSI value with the warning of `No_Combination` found in the Flag\_Warning column.

2. **Sig\_rank\_order\_log.txt**  
Log file of rank-order analysis. The program first gets the log2FC adjusts how it sorts the samples based on the +/- log2FC. The directionality of sorting is recorded in the rank log file (`Forward_Order`/`Reverse_Order`). This helps maintain overall rank position of the samples based on the directional change in expression. For example, when the FC is positive, the rank order is lowest to highest, but for negative FC the order of samples is reversed (highest to lowest). 

3. **Sig\_gene\_log.txt**  
Log file of IGSI analysis. This file is a log file of all the various combinations of the samples that get examined with corresponding IGSI score.  To calculate the IGSI, the programs determines if the lowest normalized gene count from the group with the higher average normalized gene count is higher than the highest normalized gene count in the other group. If the ranges overlap, samples are removed from either group until an IGSI is found or a cut-off threshold is breached and no IGSI is calculated.  
IGSI is calculated with the following equation:
<pre>
IGSI = (Group\_1\_Samples - Samples\_Excluded)/(Group\_1\_Samples) + 
           (Group\_2\_Samples - Samples\_Excluded)/(Group\_2\_Samples)
</pre>
For example, if no samples are removed from either grouping the score will be 2 and if 50% of samples are removed from both groups the score is 1. Once an IGSI is calculated, a new fold-change is calculated based on the samples included in the groups. The combination of samples that has the highest IGSI and highest fold-change is recorded and passed along for further analysis. During the IGSI analysis, the log file records all the combinations analyzed. For example, a user might see the following in the log file: `Group1_4+Group2_3`. This means 4 samples from Group 1 were compared to 3 samples from Group 2. If no IGSI can be calculated, the program records `No_Combination` in the log file. If a valid IGSI is calculated, the value is recorded. 

4. **Sig\_summary\_report.txt**  
Summary report of exclusion counts and rank-order analysis results for diagnosing anomalies in individual samples. The final file produced from the analysis is the summary report file. This file contains three reports: (i) a summary report showing the number of times each sample was excluded from the IGSI analysis along with its exclusion percentage (100\*(# of times samples was excluded from final gene combination)/(# of genes with passing IGSI value)\). Directly below the table are some statistics  about the dataset that help assist when analyzing the number of times a sample was excluded, (ii) a table of the average rank for each sample across all genes, and (iii) a report showing the number of times each sample appeared at each rank.

5. **FPKM\_Samples.txt**  (Only produced when using the "Cuffdiff" setting)  
Reorganized flat file of FPKM data of all genes and samples. File was created using the read\_groups.info and the genes.read\_group\_tracking file. The read\_groups.info file consists of original samples IDs with a new IDs assigned by Cuffdiff. Using this "key," the program then uses it to parse through the genes.read\_group\_tracking file, which consists of FPKM values determined by Cuffdiff, but reported in a slightly difficult manner for a user to interpret easily.

*NOTE:* Intermediate files and log files may appear redundant, but are actually important for traceability and trouble-shooting purposes. 


## References


Trapnell, C., Hendrickson, D. G., Sauvageau, M., Goff, L., Rinn, J. L., & Pachter, L. (2013). Differential analysis of gene regulation at transcript resolution with RNA-seq. Nature biotechnology, 31(1), 46-53.


## Versions of ARB
ARB2.0.0
-Enabled program to accept non-Cuffdiff datasets using template files, renamed parameter file input, removed outlier functions and removed fold-change calculations for reported IGSI groups.

ARB1.0.0 
-Fully running code


