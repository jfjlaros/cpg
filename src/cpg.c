#include <stdio.h>  // printf(), malloc() free().
#include <stdlib.h> // atoi(), atof().
#include <string.h> // strcmp().
#include <math.h>   // fmaxf() (use -lm when compiling).                    
#include <expat.h>  // XML functions (use -lexpat when compiling), provided
                    // by libexpat1-dev.

//
// Definitions.
//

#define BUFFSIZE 8192         // Buffer size for the XML parser.
#define REF_BUILD "37:GRCh37" // Use this reference sequence.


//
// Typedefs.
//

//
// Structure to hold crucial information of the candidates.
//
typedef struct _candidate {
  int pos,     // The position of the SNP.
      id;      // The dbSNP rs number.
  char nuc;    // The nucleotide of interest: C means (C, A), G means (G, T).
  float freq1, // The frequency of the first nucleotide (C or G).
        freq2; // The frequency of the second nucleotide (A or T).
} candidate;//candidate


//
// Global variables.
//

char Buff[BUFFSIZE];         // Buffer for the XML parser.
candidate *candidate_list = NULL; // List of (C, A) and (G, T) candidates.
float freqs[4],              // Current maximum frequencies for A, C, G and T.
      FreqThreshold = 0.0;   // The minimum frequency threshold.
int candidate_list_size = 0, // Size of candidate_list.
    snp_location,            // Location of the SNP under scrutiny.
    snp_id,                  // ID of the SNP under scrutiny.
    orientation;             // Variable to compensate for reverse complement.


//
// Function prototypes.
//

void add_candidate(char, float, float); // Add an entry to the candidate list.
static void XMLCALL start(void *, const char *, const char **), // XML start
                                                                // element.
                    end(void *data, const char *el);            // XML end.
                                                                // element.
int compar(const void *, const void *); // Comparison function for 
                                        // candidate_list.


//
// Functions.
//

//
// Compare two candidates based on their position.
//
// Arguments:
//   const void *a ; First candidate.
//   const void *b ; Second candidate.
//
// Returns:
//   int ; Smaller than 0 if ai.pos < pi.pos,
//         Larger than 0 if ai.pos > pi.pos,
//         0 if ai.pos = pi.pos.
//
int compar(const void *a, const void *b) {
  candidate ai = *(const candidate *)a,
            bi = *(const candidate *)b;

  return ai.pos - bi.pos;
}//compar

//
// Add an entry to the candidate list.
//
// Arguments:
//   char nuc    ; The nucleotide of interest.
//   float freq1 ; The frequency of the first nucleotide (C or G).
//         freq2 ; The frequency of the second nucleotide (A or T).
//
// Global variables:      
//   int snp_location ; The position of the SNP.
//       snp_id       ; The id of the SNP.
//
// Global variables (altered):
//   candidate_list_size ; Will be increased by one.
//   candidate_list      ; Will be expanded with a new candidate.
//
void add_candidate(char nuc, float freq1, float freq2) {
  candidate_list_size++;
  candidate_list = (candidate *)realloc(candidate_list, 
                                        candidate_list_size * 
                                        sizeof(candidate));
  candidate_list[candidate_list_size - 1].pos = snp_location;
  candidate_list[candidate_list_size - 1].id = snp_id;
  candidate_list[candidate_list_size - 1].nuc = nuc;
  candidate_list[candidate_list_size - 1].freq1 = freq1;
  candidate_list[candidate_list_size - 1].freq2 = freq2;
}//add_candidate

//
// This function is called each time a start element is encountered. Depending
// on the tag, different operations are performed:
// - SnpInfo    ; Reset the frequency table freqs[].
// - SnpLoc     ; Check if the build is REF_BUILD (defined above) and use that
//                start tag as the location of the SNP.
// - AlleleFreq ; Check if the frequency is higher than the ones previously 
//                found. Store the maximum.
//
// Arguments:
//   void *data        ; Not used.
//   const char *el    ; The element under scrutiny.
//              **attr ; Attributes of the element (if any).
//
// Global variables (altered):
//   float freqs[]    ; Altered on SnpInfo or AlleleFreq.
//   int snp_location ; Altered on SnpLoc.
//       snp_id       ; Altered on SnpInfo.
//
static void XMLCALL start(void *data, const char *el, const char **attr) {
  int i; // General iterator.

  if (!strcmp(el, "SnpInfo")) { // The SnpInfo start element.
    for (i = 0; strcmp(attr[i], "rsId"); i += 2);
    snp_id = atoi(attr[i + 1]); // Store the ID of the SNP.
    for (i = 0; i < 4; i++)     // Reset the recorded frequencies.
      freqs[i] = 0.0;
  }//if

  if (!strcmp(el, "SnpLoc")) {             // The SnpLoc start element.
    // Search for the genomicAssembly tag.
    for (i = 0; strcmp(attr[i], "genomicAssembly"); i += 2);
    if (!strcmp(attr[i + 1], REF_BUILD)) { // Is it the build we want?
      // Search for the start tag.
      for (i = 0; attr[i] && strcmp(attr[i], "start"); i += 2); 
      if (attr[i]) 
        snp_location = atoi(attr[i + 1]);    // Store the position.
      else
        snp_location = 0;
    }//if
  }//if

  if (!strcmp(el, "SsInfo")) {     // Compensate for reverse complement.
    orientation = 0;
    for (i = 0; strcmp(attr[i], "ssOrientToRs"); i += 2);
    if (!strcmp(attr[i + 1], "rev"))
      orientation = 3;
  }//if

  if (!strcmp(el, "AlleleFreq")) { // The AlleleFreq start element.
    switch (attr[1][0]) {          // Store the frequencies if applicable.
      case 'A':
        freqs[abs(orientation - 0)] = fmaxf(freqs[abs(orientation - 0)], 
                                            atof(attr[3]));
        break;
      case 'C':
        freqs[abs(orientation - 1)] = fmaxf(freqs[abs(orientation - 1)], 
                                            atof(attr[3]));
        break;
      case 'G':
        freqs[abs(orientation - 2)] = fmaxf(freqs[abs(orientation - 2)], 
                                            atof(attr[3]));
        break;
      case 'T':
        freqs[abs(orientation - 3)] = fmaxf(freqs[abs(orientation - 3)], 
                                            atof(attr[3]));
        break;
    }//switch
  }//if
}//start


//
// This function is called each time an end element is encountered.
//
// Arguments:
//   void *data        ; Not used.
//   const char *el    ; The element under scrutiny.
//
static void XMLCALL end(void *data, const char *el) {
  if (!strcmp(el, "SnpInfo")) { // The SnpLoc end element.
    if ((freqs[1] >= FreqThreshold) && (freqs[3] >= FreqThreshold))
      add_candidate('C', freqs[1], freqs[3]);
    if ((freqs[0] >= FreqThreshold) && (freqs[2] >= FreqThreshold))
      add_candidate('G', freqs[0], freqs[2]);
  }//if
}//end

//
// Entry point of the program.
//
// Arguments:
//   int argc    ; Number of arguments in argv.
//   char **argv ; Argument list:
//               ; argv[1] is an optional frequency threshold.
//
// Returns:
//   int ; 0 on no error, 1 on error.
//
int main(int argc, char **argv) {
  XML_Parser p = XML_ParserCreate(NULL); // The parser object.
  int done = 0,      // End of file indicator.
      len,           // Length of the read.
      i,             // General iterator.
      j;             // General iterator.

  if (!p) {
    fprintf(stderr, "Couldn't allocate memory for parser\n");
    return 1;
  }//if

  if (argv[1]) // The optional threshold argument.
    FreqThreshold = atof(argv[1]) / 100;

  XML_SetElementHandler(p, start, end);

  while (!done) { // Read the input and parse it.
    len = fread(Buff, 1, BUFFSIZE, stdin);
    if (ferror(stdin)) {
      fprintf(stderr, "Read error\n");
      return 1;
    }//if
    done = feof(stdin);
    if (XML_Parse(p, Buff, len, done) == XML_STATUS_ERROR) {
      fprintf(stderr, "Parse error at line %d:\n%s\n",
              (int)XML_GetCurrentLineNumber(p),
              XML_ErrorString(XML_GetErrorCode(p)));
      return 1;
    }//if
  }//while

  // Sort the candidates based on their position.
  qsort(candidate_list, candidate_list_size, sizeof(candidate), compar);

  // Scan the list for CpG islands.
  for (i = 1; i < candidate_list_size; i++) {
    for (j = i - 1; 
         j > 0 && candidate_list[i].pos == candidate_list[j].pos; j--);
    if ((candidate_list[i].nuc == 'G') && (candidate_list[j].nuc == 'C') &&
        candidate_list[i].pos == candidate_list[j].pos + 1)
      printf("%i %i %i %f %f %f %f\n", candidate_list[j].pos, 
                                       candidate_list[j].id, 
                                       candidate_list[i].id, 
                                       candidate_list[j].freq1, 
                                       candidate_list[j].freq2,
                                       candidate_list[i].freq1, 
                                       candidate_list[i].freq2);
  }//for                                         

  free(candidate_list);

  return 0;
}//main
