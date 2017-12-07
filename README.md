# ARB: Analysis of RNA-Seq Data Behavior Tool
***
*NOTE:* Program submitted to BMC Bioinformatics for review

ARB is a Python-based tool for assessment of data behavior and prioritization of differentially expressed genes, using output files from Cuffdiff, a program in the Cufflinks suite (Trapnell 2010).  

Code was developed in the Abasht Laboratory at the University of Delaware under the supervision of [Dr. Behnam Abasht](http://canr.udel.edu/faculty/behnam-abasht/).

Please contact Dr. Behnam Abasht (abasht@udel.edu) with any questions or concerns. 


## Program Overview

### Step 1. Parsing FPKM Data

The program first parses through the read\_groups.info and genes.read\_group\_tracking files produced by Cuffdiff to create a more standardized file of FPKM values. The parsed file is then filtered for statistically significant differentially expressed genes using Cuffdiff’s gene\_exp.diff file.



### Step 2. Inter-Group Separation Index (IGSI) Analysis

The parsed FPKM data is then used to analyze each gene’s overall data behavior among samples between Group 1 and Group 2 (groups are referred to as conditions in Cuffdiff) by calculating the inter-group separation index (IGSI). To calculate IGSI, the program determines if the lowest value from the group with the higher average FPKM is higher than the highest value in the other group. In other words, it determines if the two groups/conditions have overlapping ranges. If the ranges overlap, samples are removed from either group until an optimal IGSI is found or until the exclusion ceiling is breached (exclusion ceiling is set by the user). 

### Step 3. IQR Outlier Analysis

Outliers are identified using Interquartile Range (IQR) (Tukey 1977). Samples are examined on a per-gene basis to identify moderate outliers (vale = IQR x 1.5-3.0) and extreme outliers (value > IQR x 3.0) (Barbato 2010). To reduce skewness of the data, the samples are log2 transformed. Prior to transformation, a value of 1 is added to all samples to avoid taking a log of zero. 

### Step 4. Rank-Order Analysis

The overall rank-order of samples per gene is calculated and averaged across FPKM values. This helps to identify samples that consistently show higher or lower expression across genes. Additionally, the tallies of samples per rank are calculated, to show the overall spread of the data . The program automatically adjusts for +/- log2(fold change) value reported by Cuffdiff, so the order of the samples is maintained between the groups. When fold-change is positive, the order is lowest to highest. When the fold-change is negative, the order of samples is reversed highest to lowest. 


## Installation and Use

Download the program file from website and place in local working directory on computer or server.

### Required input files from Cuffdiff analysis:

* read\_groups.info (required to match sample IDs with IDs used by Cuffdiff program) 
* genes.read\_group\_tracking (supplies the FPKM values)
* gene\_exp.diff (required to filter genes reported as statistically significant)


### System requirements

**python3.6** + the following modules 

* math
* sys
* collections
* os.path
* numpy (note NumPy may need to be installed and does not come standard with python)


### To run the program locally using Python3.6 with a parameter file

Supply the following parameters in a file entitled Analysis\_CD\_Data\_Parameter\_File.txt in the same directory as the python program (just copy and paste the three lines below):


<pre>
Parameter File (space separated)

--Group_1_Name &lt;Name_of_Group_1&gt;

--Group_2_Name &lt;Name_of_Group_2&gt;
</pre>


### To run the program from a HPC environment using a submission file

<pre>
python Analysis_Cuff_Diff_Data.Version6.2.0.py \

--Group_1_Name &lt;Name_of_Group_1&gt; \

--Group_2_Name &lt;Name_of_Group_2&gt;
</pre>


### To run the program on the command line

```
python Analysis_Cuff_Diff_Data.Version6.2.0.py --Group_1_Name <Name_of_Group_1> --Group_2_Name <Name_of_Group_2> 
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

`--Group_1_Name` **REQUIRED** Name of the first group of samples analyzed by Cuffdiff aka the "first condition" (see read_groups.info file)
 
`--Group_2_Name` **REQUIRED** Name of the second group of samples analyzed by Cuffdiff aka the "condition" (see read_groups.info file)

`--Group_1_Max_Exclusion` The maximum number of samples to exclude from group 1, cannot exceed number of samples in group (default = 1)

`--Group_2_Max_Exclusion` The maximum number of samples to exclude from group 2, cannot exceed number of samples in group (default = 1)

`--Fold_Change_Cutoff` Threshold for fold-change between the two groups (default = 1.3)

`--Statistically_Significant_Data_Only` Filter the Cuffdiff dataset for only significant genes or not (default = true)

`--Input_File_Location` Location of Cuffdiff data ---REMEMBER forward slash / (default = working directory)

`--Output_File_Location` Location to write output data ---REMEMBER forward slash / (default = working directory)




### Example submission using all parameters 
*NOTE:* Assuming submission in HPC environment (aka UNIX), hence the backslashes.


<pre>
python Analysis_Cuff_Diff_Data.Version6.2.0.py \

--Group_1_Name Blue_Chickens \

--Group_2_Name Green_Chickens \

--Group_1_Max_Exclusion 1 \

--Group_2_Max_Exclusion 1 \

--Fold_Change_Cutoff 1.3 \

--Statistically_Significant_Data_Only true \

--Input_File_Location your/favorite/directory/ \

--Output_File_Location your/favorite/output/directory/ \
</pre>



## ARB Output Files 

### Example list of output files (using default settings):

* FPKM\_Samples.txt (reorganized flat file of FPKM data of all genes and samples)
* Sig\_gene\_log.txt (shows combinations of IGSI analysis)
* Sig\_IGSI\_FoldChange.txt (final results of the IGSI analysis)
* Sig\_gene\_exp\_IGSI.txt (merged results from IGSI analysis with CuffDiff gene_exp.diff data)
* Sig\_Outlier\_log.txt (log file of outlier analysis)
* Sig\_rank\_order\_log.txt (log file of ranking analysis)
* Sig\_summary\_report.txt (summary report of all analyses performed on dataset)


*NOTE:* Intermediate files and log files may appear redundant, but are actually important for traceability and trouble-shooting purposes. 

### Details of output files

#### FPKM\_Samples.txt 

Reorganized flat file of FPKM data of all genes and samples, no genes were filtered and consists of all the original data. File was created using the read\_groups.info and the genes.read\_group\_tracking file. The read\_groups.info file consists of original samples IDs with a new IDs assigned by Cuffdiff. Using this "key," the program then uses it to parse through the genes.read\_group\_tracking file, which consists of FPKM values determined by Cuffdiff, but reported in a slightly difficult manner for a user to interpret easily.  
   


#### Sig\_gene\_log.txt 

This file is a log file of all the various combinations of the samples that get examined, with corresponding IGSI score, fold-change, and verdict for analysis.  To calculate the IGSI, the programs determines if the lowest FPKM value from the group with the higher average FPKM is higher than the highest FPKM value in the other group. If the ranges overlap, samples are removed from either group until an IGSI is found or a cut-off threshold is breached and no IGSI is calculated. 

IGSI is calculated with the following equation:

![equation](http://www.sciweavers.org/tex2img.php?![equation](http://www.sciweavers.org/tex2img.php?eq=IGSI%20%3D%20%5Cfrac%7B%5Ctextup%7BGroup%201%20Samples%7D%20-%20%5Ctextup%7BSamples%20Excluded%7D%7D%7B%5Ctextup%7BGroup%201%20Samples%7D%7D%20%2B%20%5Cfrac%7B%5Ctextup%7BGroup%202%20Samples%7D-%5Ctextup%7BSamples%20Excluded%7D%7D%7B%5Ctextup%7BGroup%202%20Samples%7D%7D&bc=White&fc=Black&im=png&fs=12&ff=arev&edit=0)

*NOTE:* For example no samples removed from either grouping, the score will be 2 and if 50% of samples are removed from both groups the score is 1.

Once an IGSI is calculated, a new fold-change is calculated based on the samples included in the groups. The combination of samples that has the highest IGSI and highest fold-change is recorded and passed along for further analysis.

During the IGSI analysis, the log file records all the combinations analyzed. For example, a user might see the following in the log file: 
`Group1_4+Group2_3`. This means 4 samples from Group 1 were compared to 3 samples from Group 2. If no IGSI can be calculated, the program records `No_Combination` in the log file. If a valid IGSI is calculated, the value is recorded along with the new fold-change. If the data has no issues, it will contain the flag "Pass". Also, if the fold-change meets the fold-change threshold, the program records `Fold_Change_Sig`.

Example: `1.75|1.9094345615897184|Pass|Fold_Change_Sig`


If group 1's average FPKM is zero, the FPKM will be be adjusted to 0.001, allowing for FC to be calculated, but with the following flag `FLAG-Mean_Grp_1_Zero`. If group 2's average FPKM is zero, group 2's FPKM will be adjusted to 0.001, allowing for FC calculations, but with the following flag `FLAG-Mean_Grp_2_Zero`. For both of these situations the data will be passed onto the rest of the program for further analysis, but with these flag warnings for the user.   

If it fails the fold-change cutoff, it will instead show `Fold_Change_Not_Sig` for the last parameter`.



#### Sig\_IGSI\_FoldChange.txt

This file is the final report from the IGSI analysis for each gene. It contains information on the combination of samples that had the highest IGSI and fold-change values for each gene. The report also contains further statistical calculations based on these new groups of samples. 

The first new calculation performed in the report is the coefficient of variation (CV),  which is calculated for each "new" combination of samples. Next is the Flag\_Warning column. If  group 1 and group 2 have no issues says, this column says `Pass`. If group 1's average FPKM value was zero, the column contains the warning `FLAG-Mean_Group_1_Zero`. If group 2's average FPKM value was zero, the column displays `FLAG-Mean_Group_2_Zero`. The following columns contain standard fold-change calculations, which are self-explanatory based on their column titles. 

Finally, the report contains information on the number of samples excluded from either group 1 or 2 in the IGSI analysis, including their sample IDs.

If no IGSI value was found for that gene, it is recorded as "<" the lowest possible IGSI value with the warning of `No_Combination` found in the Flag\_Warning column.



#### Sig\_gene\_exp\_IGSI.txt 

This file consists of the merged results of the statistically significant genes from Cuffdiff gene\_exp.diff data and the IGSI final reports file (previously described). This file was created to help users sort through the information found in both files without having to manually merge the files, which greatly assists with downstream data analysis. 



#### Sig\_Outlier\_log.txt 

This file contains an extensive log of how the outliers were determined for each gene. The outlier log file also includes the IGSI calculated for each gene in case a user wants to see if there is a correlation between IGSI and outlier status.  

To calculate outliers, all samples are combined into one group and the FPKM values are corrected by adding 1 (to avoid taking log of 0) and then log2 tranforming them (to reduce skewness). These corrected values are reported in the log file. The outlier log file also contains information on the specific quartiles of the data and whisker ranges . There are two whisker ranges recorded in the file. The first whisker range refers to IQR x 1.5 and the second whisker range, labeled "extreme whisker range" refers to IQR x 3.0. Finally the outlier log records all the samples classified as an outlier, with modified FPKM value. The outliers are classified as either moderate  or extreme, based on which IQR range it exceeds. If no outliers are found, it is flagged as `No_Outlier_Data_Identified`.



#### Sig\_rank\_order\_log.txt 

This file is a log file of the rank orders for each gene. The program first gets the log2FC from Cuffdiff's gene\_exp.diff file. This log2FC information is recorded in the rank order log file. Using this information the program automatically adjusts how it sorts the samples based on the +/- log2FC. The directionality of sorting is recorded in the rank log file (`Forward_Order`/`Reverse_Order`). This helps maintain overall rank position of the samples based on the directional change in expression. For example, when the FC is positive, the rank order is lowest to highest, but for negative FC the order of samples is reversed (highest to lowest). 







#### Sig\_summary\_report.txt 

The final file produced from the analysis is the summary report file. This file contains four reports: 



1. **Summary Report of Samples Excluded from IGSI Analysis**: Reports the number of times each sample was excluded from the IGSI analysis along with its exclusion percentage (100*(# of times samples was excluded from final gene combination)/(# of genes with passing IGSI value)\). Directly below the table are some statistics  about the dataset that help assist when analyzing the number of times a sample was excluded.
2. **Outlier Summary Report**: Shows the number of times each sample was classified as a moderate or extreme outlier, outlier total and overall outlier total percentage (total\_count/no.\_of\_genes * 100).
3. **Average Rank Order Summary Report**: A table of the average rank for each sample across all genes.
4. **Tallied Summary Report of Samples Per Rank**: Shows the number of times each sample appeared at each rank.


## References

Barbato G, Barini EM, Genta G, Levi R. (2011). Features and performance of some outlier detection methods. Journal of Applied Statistics, 38, 2133–49.

Trapnell, C., Hendrickson, D. G., Sauvageau, M., Goff, L., Rinn, J. L., & Pachter, L. (2013). Differential analysis of gene regulation at transcript resolution with RNA-seq. Nature biotechnology, 31(1), 46-53.

Tukey, J. W. (1977). Exploratory data analysis. Reading, MA: Addison-Wesley Pub. Co. 
