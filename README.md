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

If you want full info about each genome id, you can also grab the woltka metadata:
https://biocore.github.io/wol/data/taxonomy/ncbi/metadata.tsv.bz2

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

6.  
