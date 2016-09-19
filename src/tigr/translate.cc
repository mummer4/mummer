#include <mummer/translate.hh>


long int Translate_DNA
(const char * A, int dnaseq_len, char * tA, int Frame)

     // function to translate dna sequence to aminoacid sequence
     // uses esttrans' headers and algo
     // A as read in from gene.h Read_String
     // tA should be an empty string malloced to atleast (Len A / 3)
     // frame is 1,2,3,4,5,6
     // returns new (strlen(tA+1)) or -1 on error

{
  int         dna_int;
  int         aaseq_len;
  int         aa_index;
  const char *dna_seq, *dna_ptr, *dna_end;
  char       *aa_ptr;

  aa_ptr =  tA + 1;
  dna_seq = A + 1;
  dna_end = A + dnaseq_len;

  if ( Frame >= 1  &&  Frame <= 3 )
    {
      dna_end -= 2;
      dna_ptr = dna_seq + Frame - 1;
      for (; dna_ptr <= dna_end;)
	{
	  dna_int = transdna [ int(*dna_ptr++) ];
	  if ( dna_int == BAD_PEP_CHAR ) {
	    fprintf(stderr,"WARNING: Forcing unrecognized DNA char to N\n");
	    dna_int = DNA_XN;
	  }
	  aa_index = dna_int;
	  aa_index <<= 4;

	  dna_int = transdna [ int(*dna_ptr++) ];
	  if ( dna_int == BAD_PEP_CHAR ) {
	    fprintf(stderr,"WARNING: Forcing unrecognized DNA char to N\n");
	    dna_int = DNA_XN;
	  }
	  aa_index += dna_int;
	  aa_index <<= 4;
	  
	  dna_int = transdna [ int(*dna_ptr++) ];
	  if ( dna_int == BAD_PEP_CHAR ) {
	    fprintf(stderr,"WARNING: Forcing unrecognized DNA char to N\n");
	    dna_int = DNA_XN;
	  }
	  aa_index += dna_int;
	  
	  *(aa_ptr++) = universal[aa_index];
	}
      
      *aa_ptr = '\0';
      aaseq_len = strlen ( tA + 1 );
    }
  else if ( Frame >= 4  &&  Frame <= 6 )
    {
      Frame -= 3;

      dna_seq += 2;
      dna_ptr = dna_end - Frame + 1;
      for (; dna_ptr >= dna_seq;)
	{
	  dna_int = compdna [ transdna [ int(*dna_ptr--)] ];
	  if ( dna_int == BAD_PEP_CHAR ) {
	    fprintf(stderr,"WARNING: Forcing unrecognized DNA char to N\n");
	    dna_int = DNA_XN;
	  }
	  aa_index = dna_int;
	  aa_index <<= 4;
	  
	  dna_int = compdna [ transdna [ int(*dna_ptr--)] ];
	  if ( dna_int == BAD_PEP_CHAR ) {
	    fprintf(stderr,"WARNING: Forcing unrecognized DNA char to N\n");
	    dna_int = DNA_XN;
	  }
	  aa_index += dna_int;
	  aa_index <<= 4;
	  
	  dna_int = compdna [ transdna [ int(*dna_ptr--)] ];
	  if ( dna_int == BAD_PEP_CHAR ) {
	    fprintf(stderr,"WARNING: Forcing unrecognized DNA char to N\n");
	    dna_int = DNA_XN;
	  }
	  aa_index += dna_int;
	  
	  *(aa_ptr++) = universal[aa_index];
	}
      
      *aa_ptr = '\0';
      aaseq_len = strlen ( tA + 1 );
    }
  else
    aaseq_len = -1;
     

  return aaseq_len;
}

