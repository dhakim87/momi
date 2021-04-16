# momi
Molecular mimicry scripts tailored to search for myelin mimics in the gut microbiome

1.  Install dependencies

You will need: Entrez Direct (edirect), either use script install-edirect here, or search for latest ncbi docs on installation (Try here: https://www.ncbi.nlm.nih.gov/books/NBK179288/)

You will need: NetMHCIIpan, Try downloading here: https://services.healthtech.dtu.dk/software.php

You will need: A conda environment to install python dependencies

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

		python download_proteomes.py

		