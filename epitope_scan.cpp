#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdexcept>

//MVAFKGVWT ISIAKSVFY VFYIAAQCL ILYLVTPPS VTMVHGNLT FTIFASCDS IALAIGFSV AINYTGASM MGNWENHWI WIYWVGPII FKEAFSKAA VVHVIDVDR VHVIDVDRG LIRLFSRDA LDVMASQKR IGRFFGGDR TAHYGSLPQ VHFFKNIVT FKNIVTPRT IFKLGGRDS IETYFSKNY INVIHAFQY YVIYGTASF FFLYGALLL LLLAEGFYT FGDYKTTIC VYIYFNTWT YIYFNTWTT VLPWNAFPG WNAFPGKVC FHLFIAAFV IAAFVGAAA IAATYNFAV IRALVGDEV VGWYRPPFS VHLYRNGKD IVPVLGPLV
const int FLANKING = 6;
const int EPITOPE_LEN = 9;
const int MAX_PROTEIN_NAME_LEN = 2048;
const int NUM_AMINO = 21;
const int BLOSUM_THRESH = 25;

const int BLOSUM62X[] = 
{
 4,  0, -2, -1, -2,  0, -2, -1, -1, -1, -1, -2, -1, -1, -1,  1,  0,  0, -3, -2, -1,
 0,  9, -3, -4, -2, -3, -3, -1, -3, -1, -1, -3, -3, -3, -3, -1, -1, -1, -2, -2, -1,
-2, -3,  6,  2, -3, -1, -1, -3, -1, -4, -3,  1, -1,  0, -2,  0, -1, -3, -4, -3, -1,
-1, -4,  2,  5, -3, -2,  0, -3,  1, -3, -2,  0, -1,  2,  0,  0, -1, -2, -3, -2, -1,
-2, -2, -3, -3,  6, -3, -1,  0, -3,  0,  0, -3, -4, -3, -3, -2, -2, -1,  1,  3, -1,
 0, -3, -1, -2, -3,  6, -2, -4, -2, -4, -3,  0, -2, -2, -2,  0, -2, -3, -2, -3, -1,
-2, -3, -1,  0, -1, -2,  8, -3, -1, -3, -2,  1, -2,  0,  0, -1, -2, -3, -2,  2, -1,
-1, -1, -3, -3,  0, -4, -3,  4, -3,  2,  1, -3, -3, -3, -3, -2, -1,  3, -3, -1, -1,
-1, -3, -1,  1, -3, -2, -1, -3,  5, -2, -1,  0, -1,  1,  2,  0, -1, -2, -3, -2, -1,
-1, -1, -4, -3,  0, -4, -3,  2, -2,  4,  2, -3, -3, -2, -2, -2, -1,  1, -2, -1, -1,
-1, -1, -3, -2,  0, -3, -2,  1, -1,  2,  5, -2, -2,  0, -1, -1, -1,  1, -1, -1, -1,
-2, -3,  1,  0, -3,  0,  1, -3,  0, -3, -2,  6, -2,  0,  0,  1,  0, -3, -4, -2, -1,
-1, -3, -1, -1, -4, -2, -2, -3, -1, -3, -2, -2,  7, -1, -2, -1, -1, -2, -4, -3, -1,
-1, -3,  0,  2, -3, -2,  0, -3,  1, -2,  0,  0, -1,  5,  1,  0, -1, -2, -2, -1, -1,
-1, -3, -2,  0, -3, -2,  0, -3,  2, -2, -1,  0, -2,  1,  5, -1, -1, -3, -3, -2, -1,
 1, -1,  0,  0, -2,  0, -1, -2,  0, -2, -1,  1, -1,  0, -1,  4,  1, -2, -3, -2, -1,
 0, -1, -1, -1, -2, -2, -2, -1, -1, -1, -1,  0, -1, -1, -1,  1,  5,  0, -2, -2, -1,
 0, -1, -3, -2, -1, -3, -3,  3, -2,  1,  1, -3, -2, -2, -3, -2,  0,  4, -3, -1, -1,
-3, -2, -4, -3,  1, -2, -2, -3, -3, -2, -1, -4, -4, -2, -3, -3, -2, -3, 11,  2, -1,
-2, -2, -3, -2,  3, -3,  2, -1, -2, -1, -1, -2, -3, -1, -2, -2, -2, -1,  2,  7, -1,
-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1
};
const int aa_map[26]= {
	0,-10000,1,2,3,4,5,6,7,-10000,8,9,10,11,-10000,12,13,14,15,16,-10000,17,18,20,19,-10000
};

//Ugh, maybe should just use asserts....
static inline int map_aa(char aa){
	int i = aa-'A';
	if (i < 0 || i > 26){
		printf("UNKNOWN AMINO: %d\n", aa);
		throw std::invalid_argument( "Unknown amino acid" );
	}
	int result = aa_map[i];
	if (result == -10000){
		printf("UNKNOWN AMINO: %c\n", aa);
		throw std::invalid_argument( "Unknown amino acid" );
	}

	return result;
}

class CircBuf
{
	public:
		int bufsize;
		char* arr;
		int index;
		int len;

		CircBuf(int bufsize)
		{
			this->bufsize = bufsize;
			this->arr = (char*)malloc(bufsize);
			this->index = 0;
			this->len = 0;
		}

		~CircBuf()
		{
			free(this->arr);
		}

		void debug_print_buf()
		{
			for (int i = 0; i < len; i++)
				printf("%c", arr[i]);
			printf("\n");
		}

		void debug_print_window()
		{
			printf("\"");
			//Start: index, End: index - 1
			if (index != len)
			{
				for (int i = index; i < bufsize; i++)
					printf("%c", arr[i]);
			}
			for (int i = 0; i < index; i++)
				printf("%c", arr[i]);
			printf("\"\n");
		}

		void copy_window(char* buffer, int offset = 0, int len = -1)
		{
			if (len < 0)
				len = bufsize;

			int start_index = 0;
			if (is_full())
				start_index = index;

			start_index = (start_index + offset) % bufsize;
			int i = start_index;
			int j = 0;
			for (;j < len; i = (i+1) % bufsize, j++)
				buffer[j] = arr[i];
			buffer[j] = '\0';
		}

		void addc(char c)
		{
			arr[index++] = c;
			if (index == bufsize)
				index = 0;
			if (len < bufsize)
				len++;
		}

		void clear()
		{
			this->index = 0;
			this->len = 0;
		}

		bool is_full()
		{
			return len == bufsize;
		}
};

char protein_name[MAX_PROTEIN_NAME_LEN];

static inline int _blosum(char c1, char c2){
	int i1 = map_aa(c1);
	int i2 = map_aa(c2);
	return BLOSUM62X[i1 + NUM_AMINO * i2];
}

int blosum(CircBuf* cb, char* epitope, int start_index, int epitope_len)
{
	int start_pos = (cb->index + start_index) % cb->bufsize;
	int end_pos = (start_pos + epitope_len) % cb->bufsize;
	int sum = 0;

	if (end_pos > start_pos)
	{
		int i = start_pos;
		int j = 0;
		for (;i < end_pos; i++,j++)
		{
			sum += _blosum(cb->arr[i], epitope[j]);
		}
	}
	else
	{
		int i = start_pos;
		int j = 0;
		for (;i < cb->bufsize; i++,j++)
		{
			sum += _blosum(cb->arr[i], epitope[j]);
		}
		for (i=0;i < end_pos; i++,j++)
		{
			sum += _blosum(cb->arr[i], epitope[j]);
		}
	}
	return sum;
}

void blosum_scan(CircBuf* cb, char** epitopes, int num_epitopes)
{
	char flanking_window[EPITOPE_LEN + 2 * FLANKING + 1];
	char window[EPITOPE_LEN + 1];


	for (int i = 0; i < num_epitopes; i++)
	{
		int blosum_score = blosum(cb, epitopes[i], FLANKING, EPITOPE_LEN);
		if (blosum_score > BLOSUM_THRESH)
		{
			cb->copy_window(flanking_window);
			cb->copy_window(window, FLANKING, EPITOPE_LEN);
			//new uuid, file_name, protein_name, flanking_window, mimic, mimicked_epitope, blosum_score
			printf("%s, %s, %s, %s, %s, %s, %d\n",
				"needUUID",
				"needfName",
				protein_name,
				flanking_window,
				window,
				epitopes[i],
				blosum_score
			);
		}
	}
}

int main(int argc, const char** argv)
{
	const char* exe_name = argv[0];
	int num_epitopes = argc-1;
	char** epitopes = new char*[num_epitopes];

	for (int i = 1; i < argc; i++)
	{
		int epi_len = strlen(argv[i]);
		if (epi_len != EPITOPE_LEN)
			throw std::invalid_argument( "Bad Epitope Length" );
		epitopes[i-1] = (char*)malloc(epi_len + 1);
		strcpy(epitopes[i-1], argv[i]);
	}

	CircBuf flanking_window(EPITOPE_LEN + 2 * FLANKING);

	const char PARSE_PROTEIN_HEADER = 0;
	const char PARSE_PROTEIN = 1;
	char state = PARSE_PROTEIN_HEADER;

	int name_index = 0;
	int c;
	int num_proteins = 0;
	while ((c = getc(stdin)) != EOF)
	{
		if (state == PARSE_PROTEIN)
		{
			if (c == '>')
			{
			    num_proteins++;
			    if (num_proteins % 1000 == 0)
			        fprintf(stderr, "%d\n", num_proteins);
				state = PARSE_PROTEIN_HEADER;
				name_index = 0;
				continue;
			}
			else if (c == '\n')
			{
				continue;
			}
			else
			{
				flanking_window.addc(c);
				if (flanking_window.is_full())
				{
					blosum_scan(&flanking_window, epitopes, num_epitopes);
				}
				continue;
			}
		}
		if (state == PARSE_PROTEIN_HEADER)
		{
			if (c == '\n')
			{
				for (int i=0; i < FLANKING; i++)
				{
					flanking_window.addc('$');
					if (flanking_window.is_full())
						blosum_scan(&flanking_window, epitopes, num_epitopes);
				}
				protein_name[name_index] = '\0';
				state = PARSE_PROTEIN;
				flanking_window.clear();
				for (int i=0; i < FLANKING; i++)
					flanking_window.addc('^');
				continue;
			}
			if (name_index == MAX_PROTEIN_NAME_LEN - 1)
				protein_name[name_index] = '\0';
			else
				protein_name[name_index++] = c;
			continue;
		}
	}

	return 0;
}
