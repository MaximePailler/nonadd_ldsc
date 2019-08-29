# nonadd_ldsc
This project is divided into two directries. *Compute_LDSC* contains all the scripts allowing to calculate the LD Scores in additive and non-additive model. *Create_Data_GWAS* contains the script in order to use these LD Scores within GWAS results.

## Compute_LDSC

### Create the database
You will find the plink script among the files. Use it to create the files you will use for the study.  
You just have to change the path in this file.  
The Python script only use .tped or .bed files. So please, do not erase these options in this script.  

### Run the Python Script
To be able to run correctly, all the .py files here must be in the same directory but you just have to run main.py  
You will be asked which type of files you want to use (.tped or .bep). Choose among these 2 options by typing '1' or '2'.

Then, you will have to write the path to the files you want to use but without typing the extension.  
For example, if you want to use _myFile.tped_, just write _/path/myFile_.

Then, you enter the minimal value of MAF you want for the SNPs inside your study.  
Please, do not enter a string or a negative value. It will return an error message and you will have to enter another value.

Finally, you will choose the size of the window for the calculation of the LD Score. The unit used here is the Morgans.  
So, if you want to see the SNPs that are at 1 centiMorgan from the SNP you are studying, just enter '0.01'.

Once the script is running, you will be able to follow its progress on the console.

The output file is a .csv one. It will be registered in the same directory as your .tped or .bed file.

### Prepare the data to use them into a GWAS
The script called *add_position.sh* will create another Id for each SNP based on its chromosome and its position in base-pair coordinates. It will be usefull to merge these results with the GWAS ones if SNPs are identified that way.

Finally, *concat_result.sh* will concatenate all your files created by this python script into one only file regrouping the data for all chromosomes.


## Create_Data_GWAS

If the results of the GWAS is divided into several files (a file for each chromosome for example), it will be necessary to regroup these files into one. That's the goal of the script called *concat_imputed.sh*.

Then, as mentioned above, you may have different formats for SNP's Id. So, in order to the merging to proceed well, you will need to split your file according to the way its identify the SNP's. For this step, you can use *split.py*.

Now the merging is possible between the GWAS files and the LD Scores ones. Use *merge.py* to do so.

Finally, you will just have to concatenate these last files into one with *concat_ld_GWAS.py*
