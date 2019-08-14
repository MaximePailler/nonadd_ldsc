# nonadd_ldsc
These are made for compute the LD Score via a certain file

## Create the database
You will find the plink script among the files. Use it to create the files you will use for the study.  
You just have to change the path in this file.  
The Python script only use .tped or .bed files. So please, do not erase these options in this script.  

## Run the Python Script
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
