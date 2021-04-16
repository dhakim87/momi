# momi
Molecular mimicry scripts tailored to search for myelin mimics in the gut microbiome

1.  Install dependencies

You will need: Entrez Direct (edirect), either use script install-edirect here, or search for latest ncbi docs on installation (Try here: https://www.ncbi.nlm.nih.gov/books/NBK179288/)

You will need: NetMHCIIpan, Try downloading here: https://services.healthtech.dtu.dk/software.php

You will need: A conda environment to install python dependencies.  You must name your conda environment "imsms" or modify the cluster job submission scripts in steps 4 and later for parallel processing.  

2.  Build list of ncbi_ids of interest.  This is everything you will download and evaluate proteome data for

File name: ncbi_taxids.txt
File Format (tab separated columns):
genome_id	ncbi_id
G000005825	398511
G000006175	456320
G000006605	306537
G000006725	160492
...

genome_id matches ids in woltka tool, ncbi_id matches ids from ncbi.  

You can download the latest list of woltka genome_ids and what they map to here:
https://biocore.github.io/wol/data/taxonomy/ncbi/

You can retrieve the taxids file, rename it, and add the column headers
https://biocore.github.io/wol/data/taxonomy/ncbi/taxids.txt

3.  Download proteome data

	3.0 - Required for large downloads: 
		Retrieve an ncbi API key, export it as an environment variable
		export NCBI_API_KEY=<whateverYourApiKeyIs>

		Entrez direct overview:
		https://www.ncbi.nlm.nih.gov/books/NBK179288/

		Retrieving an api key:
		https://ncbiinsights.ncbi.nlm.nih.gov/2017/11/02/new-api-keys-for-the-e-utilities/

	3.1 - Grab information about each ncbi_id (Relatively quick)
		python build_qinfo.py

			You can now check out your ./proteomes folder and it will be full of qinfo files containing number of proteins that need to be downloaded for each ncbi_id

	3.2 - Download proteomes (This takes minutes/days/weeks/months depending on how many proteomes you want to download and how entrez direct is feeling.  It will also silently fail at the entrez direct layer on random proteomes.)

		python download_proteome.py

	3.3 - Validate your downloads
		If you run:
		python download_proteome.py
		and everything is already downloaded correctly, it will report that it can skip all downloads and run reasonably quickly  

		Don't count on it- Entrez will likely be unable to completely download a number of the proteomes permanently, so you may have to mess with the rules in download_proteome.py to bypass certain failing downloads.  

		You can get a better understanding of what has been downloaded with:
		python validate.py

4.  Process proteome data with NetMHCIIPan (This takes a LONG TIME!)

	Note:  Processing scripts are designed for parallel use.  They parallelize work by assigning fasta files to different processes.  This is done by use of a JOB_INDEX and NUM_JOBS that is passed on the command line to each python script.  If you want to run a script in a single process, use a JOB_INDEX of 0, and a NUM_JOBS of 1.  If you want to parallelize a script, instead start X processes, using JOB_INDEX of 0,1,2,3,...,X-1, and NUM_JOBS of X.  The .sh wrapper files follow this procedure to parallelize jobs on the cluster, but may need to be modified for your cluster architecture and queuing mechanisms.  

	Note:  If you wanted to configure netMHCIIpan for alleles other than the HLA-DRB1 1501 risk factor for multiple sclerosis, you would need to modify this script.  

	Serial Usage:

	```# python <scriptname> <jobindex> <numjobs> <path_to_netMHCIIpan executable> ```
	python run_netmhciipan.py 0 1 "/Users/djhakim/netMHCIIpan/netMHCIIpan-4.0/netMHCIIpan"

	Parallel Usage:
	# FIRST: Set your path to netMHCIIpan and configure number of jobs in run_netmhciipan.sh.  
	# THEN: queue the jobs on your cluster
	qsub run_netmhciipan.sh

	For each ncbi id XXX, this will create two files
		./netmhc/XXX.mhc
		./netmhc/XXX.fin

		XXX.mhc contains the netMHCIIPan output, stripped to all weak and strong binders for each protein
		XXX.fin indicates processing successfully completed for this ncbi id, if XXX.mhc exists and XXX.fin does not exist, then either processing is still in progress, or processing failed out within netMHCIIPan

5.  Parse NetMHCIIPan files (This is relatively quick)

	Serial usage:
	python run_parse.py 0 1

	Parallel usage:
	qsub run_parse.sh

	This will parse all of the strong binders from the .mhc files and create three new files:
		./netmhc/XXX.core
		./netmhc/XXX.parseErr
		./netmhc/XXX.done

		XXX.core contains all strong binder epitopes
		XXX.parseErr contains any error messages (almost always empty)
		XXX.done indicates that the file was parsed.  If no .done exists, the script is still running or failed/was interrupted.  You should check if jobs are still queued or if XXX.parseErr has error messages in this case.  

6.  Score mimics by similarity to chosen myelin related protein

	Serial usage: Any or all of:

	python run_scoring.py 0 1 MBP
	python run_scoring.py 0 1 MOG
	python run_scoring.py 0 1 PLP1

	Parallel usage: 
	Change myelin related peptide in run_scoring.sh (pick any of MBP, MOG, PLP1)
	qsub run_scoring.sh

	This will score all of the strong binders from the .core files and create two new files:
		./netmhc/XXX.mbp
		./netmhc/XXX.mbp-scored
		or
		./netmhc/XXX.mog
		./netmhc/XXX.mog-scored
		or
		./netmhc/XXX.plp1
		./netmhc/XXX.plp1-scored

		XXX.mbp contains the blosum62 score and core epitopes from best match to worst match. 
		XXX.mbp-scored indicates that processing was successful.  Absence of this file indicates processing failed or is still in progress.  

	NOTE: To add scoring for new proteins, you must modify both run_scoring.py and score_cores.py to include epitopes for the new proteins

7.  Build mimic database

	7.1 Filling in the database requires additional metadata about each genome.  You can download the woltka metadata here:
		https://biocore.github.io/wol/data/taxonomy/ncbi/metadata.tsv.bz2
		(if link broken, check this page: https://biocore.github.io/wol/data/taxonomy/ncbi/)

		Unzip and rename the downloaded tsv file to wol_metadata.tsv

	7.2 Build the database
	
		```python make_epitope_db.py <similarity_threshold>```

		This will build ./mimics.db sqlite database containing information about all the mimics.  You need to pick a similarity threshold, this is the threshold for least similar blosum62 score that will be included in the database.  25 is a good cutoff, -99999 will keep every single epitope, of which there are millions, nearly all of which are probably useless for most purposes

Congratulations, you can now search all possible mimics in the microbiome.

You can use your favorite sqlite tools to browse the data in the database, some free ones are:
	DB Browser 
	SqliteSpy

mimics.db has three tables:
	genome
	mimic
	status

Table genome maps woltka genome_id to ncbi_id, and records genus and species according to the woltka metadata for your convenience.  

Table status summarizes current processing status.  
For every ncbi id in table genome, you can check whether or not the proteome was fully downloaded in file_status, and which, if any, of the myelin related proteins this ncbi_id was scored for in processing_status.  

If file_status is not set to "valid", then something went wrong in step 3, and the .fasta file doesn't match the .qinfo.  
If processing_status is not set to "mhc|mbp|mog|plp1", then one or more steps of the processing have failed.