*This code is developed on jupyter notebook, python 3.5

*library used are:
json
collections
itertools
numpy
pprint

*Useful files and folder
1. seq_helper.py, a helper function need to be saved it in local directory.
2. batch process, a jupyter notebook file, can be used to do the batch process.
   How to use: a.store the fastq files in the proper location(see below)
               b.fill up the 'data_label_map.txt' with proper "name label" pair(see below)
               c.Run cells in 'batch process' jupyter files, it will process all files listed 
                 in 'data_label_map.txt', and save it in "pre_process_dataset" folder.
3. single process, a jupyter notebook file, processes the desired fastq file in the proper location, 
   see comments in the code
4. This Readme file
5. raw_dataset, is a folder to hold all non processed fastq files
6. pre_process_dataset, is a fold to hold all processed txt files, for features see below

* For processing raw dataset
All fastq files need to be stored under "raw_dataset" directory.
For each fastq file stored under "raw_dataset" directory, manually assign a label to that file,
and save it inside the 'data_label_map.txt' file.
Here we are using 0,1 to distinguish the lab sources. In case there are more labs in the futurn we
can add 2,3.... Please follow the exact format in the 'data_label_map.txt', do not forget to add ',' in the middle.

* Feature stored are
1. Quality histogram, type(dictionary)
2. 4-mer histogram, type(dictionary)
3. label label, either 0 or 1, type(string) #in order to make it JSON serializable

(excluded) 4. length for each read, type(list), #not a good feature, but is useful later to see if there is error in processing, for example if some read is very short


