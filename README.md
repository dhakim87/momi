# momi
Molecular mimicry scripts tailored to search for myelin mimics in the gut microbiome

The output of the momi pipeline is mimics.db, a sqlite database of all putative molecular mimics from the gut microbiome 

1.  Install dependencies
-------------------------

You will need: 

* Entrez Direct (edirect), either use script install-edirect here, or search for latest ncbi docs on installation Try here: [https://www.ncbi.nlm.nih.gov/books/NBK179288/](https://www.ncbi.nlm.nih.gov/books/NBK179288/)

* NetMHCIIpan, Try downloading here: [https://services.healthtech.dtu.dk/software.php](https://services.healthtech.dtu.dk/software.php)

* A conda environment to install python dependencies.  You must name your conda environment "imsms" or modify the cluster job submission scripts in steps 4 and later for parallel processing.  

2.  Pick Proteomes 
-----------------------------
Build list of ncbi\_ids of interest.  This is everything you will download and evaluate proteome data for.  If you want mimics for the whole gut microbiome, your file must contain ncbi\_ids for the whole gut microbiome.  

File name: ncbi\_taxids.txt  
File Format (tab separated columns):  

	genome_id	ncbi_id  
	G000005825	398511  
	G000006175	456320  
	G000006605	306537  
	G000006725	160492  
	...

genome\_id matches ids in woltka tool, ncbi\_id matches ids from ncbi.  

You can download the latest list of woltka genome\_ids and what they map to here:
[https://biocore.github.io/wol/data/taxonomy/ncbi/](https://biocore.github.io/wol/data/taxonomy/ncbi/)

You can retrieve the taxids file, rename it, and add the column headers as indicated above
[https://biocore.github.io/wol/data/taxonomy/ncbi/taxids.txt](https://biocore.github.io/wol/data/taxonomy/ncbi/taxids.txt)

3.  Download proteome data
----------------------------------
3.0 - Required for large downloads: 

Retrieve an ncbi API key, export it as an environment variable

	export NCBI_API_KEY=<whateverYourApiKeyIs>

Entrez direct overview:
[https://www.ncbi.nlm.nih.gov/books/NBK179288/](https://www.ncbi.nlm.nih.gov/books/NBK179288/)

Retrieving an api key:
[https://ncbiinsights.ncbi.nlm.nih.gov/2017/11/02/new-api-keys-for-the-e-utilities/]
(https://ncbiinsights.ncbi.nlm.nih.gov/2017/11/02/new-api-keys-for-the-e-utilities/)

3.1 - Grab information about each ncbi_id (Relatively quick)
	
	python build_qinfo.py

You can now check out your ./proteomes folder and it will be full of qinfo files containing number of proteins that need to be downloaded for each ncbi_id

3.2 - Download proteomes (This takes minutes/days/weeks/months depending on how many proteomes you want to download and how entrez direct is feeling.  It will also silently fail at the entrez direct layer on random proteomes.)

	python download_proteome.py

3.3 - Validate your downloads: 
If you run python download_proteome.py and everything is already downloaded correctly, it will report that it can skip all downloads and run reasonably quickly  

Don't count on it- Entrez will likely be unable to completely download a number of the proteomes permanently, so you may have to mess with the rules in download_proteome.py to bypass certain failing downloads.  

You can get a better understanding of what has been downloaded with:

	python validate_downloaded_proteomes.py

4.  Process proteome data with NetMHCIIPan (This takes a LONG TIME!)
----------------------------------

Note:  Processing scripts are designed for parallel use.  They parallelize work by assigning fasta files to different processes.  This is done by use of a JOB\_INDEX and NUM\_JOBS that is passed on the command line to each python script.  If you want to run a script in a single process, use a JOB\_INDEX of 0, and a NUM\_JOBS of 1.  If you want to parallelize a script, instead start X processes, using JOB\_INDEX of 0,1,2,3,...,X-1, and NUM\_JOBS of X.  The .sh wrapper files follow this procedure to parallelize jobs on the cluster, but may need to be modified for your cluster architecture and queuing mechanisms.  

Note:  If you wanted to configure netMHCIIpan for alleles other than the HLA-DRB1 1501 risk factor for multiple sclerosis, you would need to modify  run\_netmhciipan.py.  

	Serial Usage:
	# python <scriptname> <jobindex> <numjobs> <path_to_netMHCIIpan executable> 
	python run_netmhciipan.py 0 1 "/Users/djhakim/netMHCIIpan/netMHCIIpan-4.0/netMHCIIpan"

	Parallel Usage:
	# FIRST: Edit run_netmhciipan.sh with your path to the  netMHCIIpan executable and configure number of jobs
	# THEN: queue the jobs on your cluster
	qsub run_netmhciipan.sh

For each ncbi id XXX, this will create two files:

  * ./netmhc/XXX.mhc - contains the netMHCIIPan output, stripped to all weak and strong binders for each protein
  * ./netmhc/XXX.fin - indicates processing successfully completed for this ncbi id, if XXX.mhc exists and XXX.fin does not exist, then either processing is still in progress, or processing failed out within netMHCIIPan

5.  Parse NetMHCIIPan files (This is relatively quick)
-------------------------------------
	Serial usage:
	python run_parse.py 0 1

	Parallel usage:
	qsub run_parse.sh

This will parse all of the strong binders from the .mhc files and create three new files:

* ./netmhc/XXX.core - contains all strong binder epitopes
* ./netmhc/XXX.parseErr - contains any error messages (almost always empty)
* ./netmhc/XXX.done - indicates that the file was parsed.  If no .done exists, the script is still running or failed/was interrupted.  You should check if jobs are still queued or if XXX.parseErr has error messages in this case.

6.  Score mimics by similarity to chosen myelin related protein
----------------------------------
	Serial usage: Any or all of:
	python run_scoring.py 0 1 MBP
	python run_scoring.py 0 1 MOG
	python run_scoring.py 0 1 PLP1

	Parallel usage: 
	Change myelin related peptide in run_scoring.sh (pick any of MBP, MOG, PLP1)
	qsub run_scoring.sh

This will score all of the strong binders from the .core files and create two new files:

*		./netmhc/XXX.(mbp/mog/plp1) - contains the blosum62 score and core epitopes sorted from best match to worst match
*		./netmhc/XXX.(mbp/mog/plp1)-scored - presence indicates processing was successful.  Absence indicates processing failed or is still in progress

NOTE: To add scoring for new proteins, you must modify both run_scoring.py and score_cores.py to include epitopes for the new proteins

7.  Build mimic database
---------------------------------

7.1 Filling in the database requires additional metadata about each genome.  You can download the woltka metadata here:  
[https://biocore.github.io/wol/data/taxonomy/ncbi/metadata.tsv.bz2](https://biocore.github.io/wol/data/taxonomy/ncbi/metadata.tsv.bz2)  
(if link broken, check this page: [https://biocore.github.io/wol/data/taxonomy/ncbi/](https://biocore.github.io/wol/data/taxonomy/ncbi/))

	# Unzip and rename the downloaded tsv file to wol_metadata.tsv
	mv metadata.tsv wol_metadata.tsv

7.2 Build the database
	
	python make_epitope_db.py <similarity_threshold>

This will build ./mimics.db sqlite database containing information about all the mimics.  You need to pick a similarity threshold, this is the threshold for least similar blosum62 score that will be included in the database.  25 is a good cutoff, -99999 will keep every single epitope, of which there are millions, nearly all of which are probably useless for most purposes

---

Congratulations, you can now search all possible mimics in the microbiome.

---

8. Browsing The Database
--------------------------
You can use your favorite sqlite tools to browse the data in the database, some free ones are:

* DB Browser 
* SqliteSpy

mimics.db has three tables:
	
* genome
* status
* mimic

Table genome maps woltka genome\_id to ncbi\_id, and records genus and species according to the woltka metadata for your convenience.  

Table status summarizes current processing status.  
For every ncbi id in table genome, you can check whether or not the proteome was fully downloaded in file\_status, and which, if any, of the myelin related proteins this ncbi\_id was scored for in processing\_status.  

If file\_status is not set to "valid", then something went wrong in step 3, and the .fasta file doesn't match the .qinfo.  
If processing_status is not set to "mhc|mbp|mog|plp1", then one or more steps of the processing have failed.

Table mimic is the set of all potential molecular mimics.  Each row shows the ncbi\_id and genome\_id of the species' reference proteome from which the mimic was sourced.  epitope denotes the the 9 amino acid core peptide of the mimic, mimicked shows which myelin related peptide was mimicked, and blosum shows the blosum62 similarity score of that mimic against the expected core epitope (or epitopes) for that myelin related protein.  

Since the range of blosum scores varies by protein (and theoretically the mimicked epitope as well, but that isn't recorded yet), it is recommended to order this database first by mimicked then by blosum score in descending order.  

A useful sql query might be:

	-- Find all mog mimics from best to worst match
	SELECT * FROM mimic LEFT JOIN genome USING (ncbi_id) WHERE mimicked='MOG' ORDER BY blosum DESC