//-- Values and function for DNA -> Amino Acid translation

#ifndef __TRANSLATE_HH
#define __TRANSLATE_HH

#include <cstdio>
#include <cstring>

long int Translate_DNA
(const char * A, int dnaseq_len, char * tA, int Frame);

     // function to translate dna sequence to aminoacid sequence
     // uses esttrans' headers and algo
     // A as read in from gene.h Read_String
     // tA should be an empty string, malloced to atleast (len A / 3)
     // frame is 1,2,3,4,5,6
     // returns new (strlen(A+1)) or -1 on error
     // See end of file for source.

inline long int Translate_DNA
(const char * A, char * tA, int Frame) {
  return Translate_DNA(A, (int)strlen(A + 1), tA, Frame);
}


#define BAD_PEP_CHAR	-1
#define SKIP_PEP_CHAR	-2
#define NUM_DNA_SYMBOLS	16
#define DNA_A		 0
#define DNA_C		 1
#define DNA_G		 2
#define DNA_TU		 3
#define DNA_M		 4
#define DNA_R		 5
#define DNA_W		 6
#define DNA_S		 7
#define DNA_Y		 8
#define DNA_K		 9
#define DNA_V		10
#define DNA_H		11
#define DNA_D		12
#define DNA_B		13
#define DNA_XN		14
#define DNA_dot		15

#define AA_LEN		8192
#define DNA_LEN		32768
#define AA_LINE_LEN	60

const int	compdna[NUM_DNA_SYMBOLS] = {
		DNA_TU,
		DNA_G,
		DNA_C,
		DNA_A,
		DNA_K,
		DNA_Y,
		DNA_W,
		DNA_S,
		DNA_R,
		DNA_M,
		DNA_B,
		DNA_D,
		DNA_H,
		DNA_V,
		DNA_XN,
		DNA_dot
};

const int	transdna[256] = {
		BAD_PEP_CHAR,	/*  0 NUL*/
		BAD_PEP_CHAR,	/*  1 SOH*/
		BAD_PEP_CHAR,	/*  2 STX*/
		BAD_PEP_CHAR,	/*  3 ETX*/
		BAD_PEP_CHAR,	/*  4 EOT*/
		BAD_PEP_CHAR,	/*  5 ENQ*/
		BAD_PEP_CHAR,	/*  6 ACK*/
		BAD_PEP_CHAR,	/*  7 BEL*/
		BAD_PEP_CHAR,	/*  8 BS */
		SKIP_PEP_CHAR,	/*  9 HT */
		SKIP_PEP_CHAR,	/* 10 NL */
		SKIP_PEP_CHAR,	/* 11 VT */
		SKIP_PEP_CHAR,	/* 12 NP */
		SKIP_PEP_CHAR,	/* 13 CR */
		BAD_PEP_CHAR,	/* 14 SO */
		BAD_PEP_CHAR,	/* 15 SI */
		BAD_PEP_CHAR,	/* 16 DLE*/
		BAD_PEP_CHAR,	/* 17 DC1*/
		BAD_PEP_CHAR,	/* 18 DC2*/
		BAD_PEP_CHAR,	/* 19 DC3*/
		BAD_PEP_CHAR,	/* 20 DC4*/
		BAD_PEP_CHAR,	/* 21 NAK*/
		BAD_PEP_CHAR,	/* 22 SYN*/
		BAD_PEP_CHAR,	/* 23 ETB*/
		BAD_PEP_CHAR,	/* 24 CAN*/
		BAD_PEP_CHAR,	/* 25 EM */
		BAD_PEP_CHAR,	/* 26 SUB*/
		BAD_PEP_CHAR,	/* 27 ESC*/
		BAD_PEP_CHAR,	/* 28 FS */
		BAD_PEP_CHAR,	/* 29 GS */
		BAD_PEP_CHAR,	/* 30 RS */
		BAD_PEP_CHAR,	/* 31 US */
		SKIP_PEP_CHAR,	/* 32 SP */
		BAD_PEP_CHAR,	/* 33  ! */
		BAD_PEP_CHAR,	/* 34  " */
		BAD_PEP_CHAR,	/* 35  # */
		BAD_PEP_CHAR,	/* 36  $ */
		BAD_PEP_CHAR,	/* 37  % */
		BAD_PEP_CHAR,	/* 38  & */
		BAD_PEP_CHAR,	/* 39  ' */
		BAD_PEP_CHAR,	/* 40  ( */
		BAD_PEP_CHAR,	/* 41  ) */
		BAD_PEP_CHAR,	/* 42  * */
		BAD_PEP_CHAR,	/* 43  + */
		BAD_PEP_CHAR,	/* 44  , */
		BAD_PEP_CHAR,	/* 45  - */
		DNA_dot,	/* 46  . */
		BAD_PEP_CHAR,	/* 47  / */
		BAD_PEP_CHAR,	/* 48  0 */
		BAD_PEP_CHAR,	/* 49  1 */
		BAD_PEP_CHAR,	/* 50  2 */
		BAD_PEP_CHAR,	/* 51  3 */
		BAD_PEP_CHAR,	/* 52  4 */
		BAD_PEP_CHAR,	/* 53  5 */
		BAD_PEP_CHAR,	/* 54  6 */
		BAD_PEP_CHAR,	/* 55  7 */
		BAD_PEP_CHAR,	/* 56  8 */
		BAD_PEP_CHAR,	/* 57  9 */
		BAD_PEP_CHAR,	/* 58  : */
		BAD_PEP_CHAR,	/* 59  ; */
		BAD_PEP_CHAR,	/* 60  < */
		BAD_PEP_CHAR,	/* 61  = */
		BAD_PEP_CHAR,	/* 62  > */
		BAD_PEP_CHAR,	/* 63  ? */
		BAD_PEP_CHAR,	/* 64  @ */
		DNA_A,		/* 65  A */
		DNA_B,		/* 66  B */
		DNA_C,		/* 67  C */
		DNA_D,		/* 68  D */
		BAD_PEP_CHAR,	/* 69  E */
		BAD_PEP_CHAR,	/* 70  F */
		DNA_G,		/* 71  G */
		DNA_H,		/* 72  H */
		BAD_PEP_CHAR,	/* 73  I */
		BAD_PEP_CHAR,	/* 74  J */
		DNA_K,		/* 75  K */
		BAD_PEP_CHAR,	/* 76  L */
		DNA_M,		/* 77  M */
		DNA_XN,		/* 78  N */
		BAD_PEP_CHAR,	/* 79  O */
		BAD_PEP_CHAR,	/* 80  P */
		BAD_PEP_CHAR,	/* 81  Q */
		DNA_R,		/* 82  R */
		DNA_S,		/* 83  S */
		DNA_TU,		/* 84  T */
		DNA_TU,		/* 85  U */
		DNA_V,		/* 86  V */
		DNA_W,		/* 87  W */
		DNA_XN,		/* 88  X */
		DNA_Y,		/* 89  Y */
		BAD_PEP_CHAR,	/* 90  Z */
		BAD_PEP_CHAR,	/* 91  [ */
		BAD_PEP_CHAR,	/* 92  \ */
		BAD_PEP_CHAR,	/* 93  ] */
		BAD_PEP_CHAR,	/* 94  ^ */
		BAD_PEP_CHAR,	/* 95  _ */
		BAD_PEP_CHAR,	/* 96  ` */
		DNA_A,		/* 97  a */
		DNA_B,		/* 98  b */
		DNA_C,		/* 99  c */
		DNA_D,		/*100  d */
		BAD_PEP_CHAR,	/*101  e */
		BAD_PEP_CHAR,	/*102  f */
		DNA_G,		/*103  g */
		DNA_H,		/*104  h */
		BAD_PEP_CHAR,	/*105  i */
		BAD_PEP_CHAR,	/*106  j */
		DNA_K,		/*107  k */
		BAD_PEP_CHAR,	/*108  l */
		DNA_M,		/*109  m */
		DNA_XN,		/*110  n */
		BAD_PEP_CHAR,	/*111  o */
		BAD_PEP_CHAR,	/*112  p */
		BAD_PEP_CHAR,	/*113  q */
		DNA_R,		/*114  r */
		DNA_S,		/*115  s */
		DNA_TU,		/*116  t */
		DNA_TU,		/*117  u */
		DNA_V,		/*118  v */
		DNA_W,		/*119  w */
		DNA_XN,		/*120  x */
		DNA_Y,		/*121  y */
		BAD_PEP_CHAR,	/*122  z */
		BAD_PEP_CHAR,	/*123  { */
		BAD_PEP_CHAR,	/*124  | */
		BAD_PEP_CHAR,	/*125  } */
		BAD_PEP_CHAR,	/*126  ~ */
		BAD_PEP_CHAR,	/*127 DEL*/
		BAD_PEP_CHAR,
		BAD_PEP_CHAR,
		BAD_PEP_CHAR,
		BAD_PEP_CHAR,
		BAD_PEP_CHAR,
		BAD_PEP_CHAR,
		BAD_PEP_CHAR,
		BAD_PEP_CHAR,
		BAD_PEP_CHAR,
		BAD_PEP_CHAR,
		BAD_PEP_CHAR,
		BAD_PEP_CHAR,
		BAD_PEP_CHAR,
		BAD_PEP_CHAR,
		BAD_PEP_CHAR,
		BAD_PEP_CHAR,
		BAD_PEP_CHAR,
		BAD_PEP_CHAR,
		BAD_PEP_CHAR,
		BAD_PEP_CHAR,
		BAD_PEP_CHAR,
		BAD_PEP_CHAR,
		BAD_PEP_CHAR,
		BAD_PEP_CHAR,
		BAD_PEP_CHAR,
		BAD_PEP_CHAR,
		BAD_PEP_CHAR,
		BAD_PEP_CHAR,
		BAD_PEP_CHAR,
		BAD_PEP_CHAR,
		BAD_PEP_CHAR,
		BAD_PEP_CHAR,
		BAD_PEP_CHAR,
		BAD_PEP_CHAR,
		BAD_PEP_CHAR,
		BAD_PEP_CHAR,
		BAD_PEP_CHAR,
		BAD_PEP_CHAR,
		BAD_PEP_CHAR,
		BAD_PEP_CHAR,
		BAD_PEP_CHAR,
		BAD_PEP_CHAR,
		BAD_PEP_CHAR,
		BAD_PEP_CHAR,
		BAD_PEP_CHAR,
		BAD_PEP_CHAR,
		BAD_PEP_CHAR,
		BAD_PEP_CHAR,
		BAD_PEP_CHAR,
		BAD_PEP_CHAR,
		BAD_PEP_CHAR,
		BAD_PEP_CHAR,
		BAD_PEP_CHAR,
		BAD_PEP_CHAR,
		BAD_PEP_CHAR,
		BAD_PEP_CHAR,
		BAD_PEP_CHAR,
		BAD_PEP_CHAR,
		BAD_PEP_CHAR,
		BAD_PEP_CHAR,
		BAD_PEP_CHAR,
		BAD_PEP_CHAR,
		BAD_PEP_CHAR,
		BAD_PEP_CHAR,
		BAD_PEP_CHAR,
		BAD_PEP_CHAR,
		BAD_PEP_CHAR,
		BAD_PEP_CHAR,
		BAD_PEP_CHAR,
		BAD_PEP_CHAR,
		BAD_PEP_CHAR,
		BAD_PEP_CHAR,
		BAD_PEP_CHAR,
		BAD_PEP_CHAR,
		BAD_PEP_CHAR,
		BAD_PEP_CHAR,
		BAD_PEP_CHAR,
		BAD_PEP_CHAR,
		BAD_PEP_CHAR,
		BAD_PEP_CHAR,
		BAD_PEP_CHAR,
		BAD_PEP_CHAR,
		BAD_PEP_CHAR,
		BAD_PEP_CHAR,
		BAD_PEP_CHAR,
		BAD_PEP_CHAR,
		BAD_PEP_CHAR,
		BAD_PEP_CHAR,
		BAD_PEP_CHAR,
		BAD_PEP_CHAR,
		BAD_PEP_CHAR,
		BAD_PEP_CHAR,
		BAD_PEP_CHAR,
		BAD_PEP_CHAR,
		BAD_PEP_CHAR,
		BAD_PEP_CHAR,
		BAD_PEP_CHAR,
		BAD_PEP_CHAR,
		BAD_PEP_CHAR,
		BAD_PEP_CHAR,
		BAD_PEP_CHAR,
		BAD_PEP_CHAR,
		BAD_PEP_CHAR,
		BAD_PEP_CHAR,
		BAD_PEP_CHAR,
		BAD_PEP_CHAR,
		BAD_PEP_CHAR,
		BAD_PEP_CHAR,
		BAD_PEP_CHAR,
		BAD_PEP_CHAR,
		BAD_PEP_CHAR,
		BAD_PEP_CHAR,
		BAD_PEP_CHAR,
		BAD_PEP_CHAR,
		BAD_PEP_CHAR,
		BAD_PEP_CHAR,
		BAD_PEP_CHAR,
		BAD_PEP_CHAR,
		BAD_PEP_CHAR,
		BAD_PEP_CHAR,
		BAD_PEP_CHAR,
		BAD_PEP_CHAR,
		BAD_PEP_CHAR,
		BAD_PEP_CHAR,
		BAD_PEP_CHAR,
		BAD_PEP_CHAR,
		BAD_PEP_CHAR,
		BAD_PEP_CHAR
};

const char	universal[NUM_DNA_SYMBOLS * NUM_DNA_SYMBOLS * NUM_DNA_SYMBOLS] = {
		'K','N','K','N','X','K','X','X','N','X','X','X','X','X','X','X',
		'T','T','T','T','T','T','T','T','T','T','T','T','T','T','T','X',
		'R','S','R','S','X','R','X','X','S','X','X','X','X','X','X','X',
		'I','I','M','I','I','X','I','X','I','X','X','I','X','X','X','X',
		'X','X','X','X','X','X','X','X','X','X','X','X','X','X','X','X',
		'X','X','X','X','X','X','X','X','X','X','X','X','X','X','X','X',
		'X','X','X','X','X','X','X','X','X','X','X','X','X','X','X','X',
		'X','X','X','X','X','X','X','X','X','X','X','X','X','X','X','X',
		'X','X','X','X','X','X','X','X','X','X','X','X','X','X','X','X',
		'X','X','X','X','X','X','X','X','X','X','X','X','X','X','X','X',
		'X','X','X','X','X','X','X','X','X','X','X','X','X','X','X','X',
		'X','X','X','X','X','X','X','X','X','X','X','X','X','X','X','X',
		'X','X','X','X','X','X','X','X','X','X','X','X','X','X','X','X',
		'X','X','X','X','X','X','X','X','X','X','X','X','X','X','X','X',
		'X','X','X','X','X','X','X','X','X','X','X','X','X','X','X','X',
		'X','X','X','X','X','X','X','X','X','X','X','X','X','X','X','X',
		'Q','H','Q','H','X','Q','X','X','H','X','X','X','X','X','X','X',
		'P','P','P','P','P','P','P','P','P','P','P','P','P','P','P','X',
		'R','R','R','R','R','R','R','R','R','R','R','R','R','R','R','X',
		'L','L','L','L','L','L','L','L','L','L','L','L','L','L','L','X',
		'X','X','X','X','X','X','X','X','X','X','X','X','X','X','X','X',
		'X','X','X','X','X','X','X','X','X','X','X','X','X','X','X','X',
		'X','X','X','X','X','X','X','X','X','X','X','X','X','X','X','X',
		'X','X','X','X','X','X','X','X','X','X','X','X','X','X','X','X',
		'X','X','X','X','X','X','X','X','X','X','X','X','X','X','X','X',
		'X','X','X','X','X','X','X','X','X','X','X','X','X','X','X','X',
		'X','X','X','X','X','X','X','X','X','X','X','X','X','X','X','X',
		'X','X','X','X','X','X','X','X','X','X','X','X','X','X','X','X',
		'X','X','X','X','X','X','X','X','X','X','X','X','X','X','X','X',
		'X','X','X','X','X','X','X','X','X','X','X','X','X','X','X','X',
		'X','X','X','X','X','X','X','X','X','X','X','X','X','X','X','X',
		'X','X','X','X','X','X','X','X','X','X','X','X','X','X','X','X',
		'E','D','E','D','X','E','X','X','D','X','X','X','X','X','X','X',
		'A','A','A','A','A','A','A','A','A','A','A','A','A','A','A','X',
		'G','G','G','G','G','G','G','G','G','G','G','G','G','G','G','X',
		'V','V','V','V','V','V','V','V','V','V','V','V','V','V','V','X',
		'X','X','X','X','X','X','X','X','X','X','X','X','X','X','X','X',
		'X','X','X','X','X','X','X','X','X','X','X','X','X','X','X','X',
		'X','X','X','X','X','X','X','X','X','X','X','X','X','X','X','X',
		'X','X','X','X','X','X','X','X','X','X','X','X','X','X','X','X',
		'X','X','X','X','X','X','X','X','X','X','X','X','X','X','X','X',
		'X','X','X','X','X','X','X','X','X','X','X','X','X','X','X','X',
		'X','X','X','X','X','X','X','X','X','X','X','X','X','X','X','X',
		'X','X','X','X','X','X','X','X','X','X','X','X','X','X','X','X',
		'X','X','X','X','X','X','X','X','X','X','X','X','X','X','X','X',
		'X','X','X','X','X','X','X','X','X','X','X','X','X','X','X','X',
		'X','X','X','X','X','X','X','X','X','X','X','X','X','X','X','X',
		'X','X','X','X','X','X','X','X','X','X','X','X','X','X','X','X',
		'*','Y','*','Y','X','*','X','X','Y','X','X','X','X','X','X','X',
		'S','S','S','S','S','S','S','S','S','S','S','S','S','S','S','X',
		'*','C','W','C','X','X','X','X','C','X','X','X','X','X','X','X',
		'L','F','L','F','X','L','X','X','F','X','X','X','X','X','X','X',
		'X','X','X','X','X','X','X','X','X','X','X','X','X','X','X','X',
		'*','X','X','X','X','X','X','X','X','X','X','X','X','X','X','X',
		'X','X','X','X','X','X','X','X','X','X','X','X','X','X','X','X',
		'X','X','X','X','X','X','X','X','X','X','X','X','X','X','X','X',
		'X','X','X','X','X','X','X','X','X','X','X','X','X','X','X','X',
		'X','X','X','X','X','X','X','X','X','X','X','X','X','X','X','X',
		'X','X','X','X','X','X','X','X','X','X','X','X','X','X','X','X',
		'X','X','X','X','X','X','X','X','X','X','X','X','X','X','X','X',
		'X','X','X','X','X','X','X','X','X','X','X','X','X','X','X','X',
		'X','X','X','X','X','X','X','X','X','X','X','X','X','X','X','X',
		'X','X','X','X','X','X','X','X','X','X','X','X','X','X','X','X',
		'X','X','X','X','X','X','X','X','X','X','X','X','X','X','X','X',
		'X','X','X','X','X','X','X','X','X','X','X','X','X','X','X','X',
		'X','X','X','X','X','X','X','X','X','X','X','X','X','X','X','X',
		'R','X','R','X','X','R','X','X','X','X','X','X','X','X','X','X',
		'X','X','X','X','X','X','X','X','X','X','X','X','X','X','X','X',
		'X','X','X','X','X','X','X','X','X','X','X','X','X','X','X','X',
		'X','X','X','X','X','X','X','X','X','X','X','X','X','X','X','X',
		'X','X','X','X','X','X','X','X','X','X','X','X','X','X','X','X',
		'X','X','X','X','X','X','X','X','X','X','X','X','X','X','X','X',
		'X','X','X','X','X','X','X','X','X','X','X','X','X','X','X','X',
		'X','X','X','X','X','X','X','X','X','X','X','X','X','X','X','X',
		'X','X','X','X','X','X','X','X','X','X','X','X','X','X','X','X',
		'X','X','X','X','X','X','X','X','X','X','X','X','X','X','X','X',
		'X','X','X','X','X','X','X','X','X','X','X','X','X','X','X','X',
		'X','X','X','X','X','X','X','X','X','X','X','X','X','X','X','X',
		'X','X','X','X','X','X','X','X','X','X','X','X','X','X','X','X',
		'X','X','X','X','X','X','X','X','X','X','X','X','X','X','X','X',
		'X','B','X','B','X','X','X','X','B','X','X','X','X','X','X','X',
		'X','X','X','X','X','X','X','X','X','X','X','X','X','X','X','X',
		'X','X','X','X','X','X','X','X','X','X','X','X','X','X','X','X',
		'X','X','X','X','X','X','X','X','X','X','X','X','X','X','X','X',
		'X','X','X','X','X','X','X','X','X','X','X','X','X','X','X','X',
		'X','X','X','X','X','X','X','X','X','X','X','X','X','X','X','X',
		'X','X','X','X','X','X','X','X','X','X','X','X','X','X','X','X',
		'X','X','X','X','X','X','X','X','X','X','X','X','X','X','X','X',
		'X','X','X','X','X','X','X','X','X','X','X','X','X','X','X','X',
		'X','X','X','X','X','X','X','X','X','X','X','X','X','X','X','X',
		'X','X','X','X','X','X','X','X','X','X','X','X','X','X','X','X',
		'X','X','X','X','X','X','X','X','X','X','X','X','X','X','X','X',
		'X','X','X','X','X','X','X','X','X','X','X','X','X','X','X','X',
		'X','X','X','X','X','X','X','X','X','X','X','X','X','X','X','X',
		'X','X','X','X','X','X','X','X','X','X','X','X','X','X','X','X',
		'X','X','X','X','X','X','X','X','X','X','X','X','X','X','X','X',
		'X','X','X','X','X','X','X','X','X','X','X','X','X','X','X','X',
		'X','X','X','X','X','X','X','X','X','X','X','X','X','X','X','X',
		'X','X','X','X','X','X','X','X','X','X','X','X','X','X','X','X',
		'X','X','X','X','X','X','X','X','X','X','X','X','X','X','X','X',
		'X','X','X','X','X','X','X','X','X','X','X','X','X','X','X','X',
		'X','X','X','X','X','X','X','X','X','X','X','X','X','X','X','X',
		'X','X','X','X','X','X','X','X','X','X','X','X','X','X','X','X',
		'X','X','X','X','X','X','X','X','X','X','X','X','X','X','X','X',
		'X','X','X','X','X','X','X','X','X','X','X','X','X','X','X','X',
		'X','X','X','X','X','X','X','X','X','X','X','X','X','X','X','X',
		'X','X','X','X','X','X','X','X','X','X','X','X','X','X','X','X',
		'X','X','X','X','X','X','X','X','X','X','X','X','X','X','X','X',
		'X','X','X','X','X','X','X','X','X','X','X','X','X','X','X','X',
		'X','X','X','X','X','X','X','X','X','X','X','X','X','X','X','X',
		'X','X','X','X','X','X','X','X','X','X','X','X','X','X','X','X',
		'X','X','X','X','X','X','X','X','X','X','X','X','X','X','X','X',
		'Z','X','Z','X','X','Z','X','X','X','X','X','X','X','X','X','X',
		'X','X','X','X','X','X','X','X','X','X','X','X','X','X','X','X',
		'X','X','X','X','X','X','X','X','X','X','X','X','X','X','X','X',
		'X','X','X','X','X','X','X','X','X','X','X','X','X','X','X','X',
		'X','X','X','X','X','X','X','X','X','X','X','X','X','X','X','X',
		'X','X','X','X','X','X','X','X','X','X','X','X','X','X','X','X',
		'X','X','X','X','X','X','X','X','X','X','X','X','X','X','X','X',
		'X','X','X','X','X','X','X','X','X','X','X','X','X','X','X','X',
		'X','X','X','X','X','X','X','X','X','X','X','X','X','X','X','X',
		'X','X','X','X','X','X','X','X','X','X','X','X','X','X','X','X',
		'X','X','X','X','X','X','X','X','X','X','X','X','X','X','X','X',
		'X','X','X','X','X','X','X','X','X','X','X','X','X','X','X','X',
		'X','X','X','X','X','X','X','X','X','X','X','X','X','X','X','X',
		'X','X','X','X','X','X','X','X','X','X','X','X','X','X','X','X',
		'X','X','X','X','X','X','X','X','X','X','X','X','X','X','X','X',
		'X','X','X','X','X','X','X','X','X','X','X','X','X','X','X','X',
		'X','X','X','X','X','X','X','X','X','X','X','X','X','X','X','X',
		'X','X','X','X','X','X','X','X','X','X','X','X','X','X','X','X',
		'X','X','X','X','X','X','X','X','X','X','X','X','X','X','X','X',
		'L','X','L','X','X','L','X','X','X','X','X','X','X','X','X','X',
		'X','X','X','X','X','X','X','X','X','X','X','X','X','X','X','X',
		'X','X','X','X','X','X','X','X','X','X','X','X','X','X','X','X',
		'X','X','X','X','X','X','X','X','X','X','X','X','X','X','X','X',
		'X','X','X','X','X','X','X','X','X','X','X','X','X','X','X','X',
		'X','X','X','X','X','X','X','X','X','X','X','X','X','X','X','X',
		'X','X','X','X','X','X','X','X','X','X','X','X','X','X','X','X',
		'X','X','X','X','X','X','X','X','X','X','X','X','X','X','X','X',
		'X','X','X','X','X','X','X','X','X','X','X','X','X','X','X','X',
		'X','X','X','X','X','X','X','X','X','X','X','X','X','X','X','X',
		'X','X','X','X','X','X','X','X','X','X','X','X','X','X','X','X',
		'X','X','X','X','X','X','X','X','X','X','X','X','X','X','X','X',
		'X','X','X','X','X','X','X','X','X','X','X','X','X','X','X','X',
		'X','X','X','X','X','X','X','X','X','X','X','X','X','X','X','X',
		'X','X','X','X','X','X','X','X','X','X','X','X','X','X','X','X',
		'X','X','X','X','X','X','X','X','X','X','X','X','X','X','X','X',
		'X','X','X','X','X','X','X','X','X','X','X','X','X','X','X','X',
		'X','X','X','X','X','X','X','X','X','X','X','X','X','X','X','X',
		'X','X','X','X','X','X','X','X','X','X','X','X','X','X','X','X',
		'X','X','X','X','X','X','X','X','X','X','X','X','X','X','X','X',
		'X','X','X','X','X','X','X','X','X','X','X','X','X','X','X','X',
		'X','X','X','X','X','X','X','X','X','X','X','X','X','X','X','X',
		'X','X','X','X','X','X','X','X','X','X','X','X','X','X','X','X',
		'X','X','X','X','X','X','X','X','X','X','X','X','X','X','X','X',
		'X','X','X','X','X','X','X','X','X','X','X','X','X','X','X','X',
		'X','X','X','X','X','X','X','X','X','X','X','X','X','X','X','X',
		'X','X','X','X','X','X','X','X','X','X','X','X','X','X','X','X',
		'X','X','X','X','X','X','X','X','X','X','X','X','X','X','X','X',
		'X','X','X','X','X','X','X','X','X','X','X','X','X','X','X','X',
		'X','X','X','X','X','X','X','X','X','X','X','X','X','X','X','X',
		'X','X','X','X','X','X','X','X','X','X','X','X','X','X','X','X',
		'X','X','X','X','X','X','X','X','X','X','X','X','X','X','X','X',
		'X','X','X','X','X','X','X','X','X','X','X','X','X','X','X','X',
		'X','X','X','X','X','X','X','X','X','X','X','X','X','X','X','X',
		'X','X','X','X','X','X','X','X','X','X','X','X','X','X','X','X',
		'X','X','X','X','X','X','X','X','X','X','X','X','X','X','X','X',
		'X','X','X','X','X','X','X','X','X','X','X','X','X','X','X','X',
		'X','X','X','X','X','X','X','X','X','X','X','X','X','X','X','X',
		'X','X','X','X','X','X','X','X','X','X','X','X','X','X','X','X',
		'X','X','X','X','X','X','X','X','X','X','X','X','X','X','X','X',
		'X','X','X','X','X','X','X','X','X','X','X','X','X','X','X','X',
		'X','X','X','X','X','X','X','X','X','X','X','X','X','X','X','X',
		'X','X','X','X','X','X','X','X','X','X','X','X','X','X','X','X',
		'X','X','X','X','X','X','X','X','X','X','X','X','X','X','X','X',
		'X','X','X','X','X','X','X','X','X','X','X','X','X','X','X','X',
		'X','X','X','X','X','X','X','X','X','X','X','X','X','X','X','X',
		'X','X','X','X','X','X','X','X','X','X','X','X','X','X','X','X',
		'X','X','X','X','X','X','X','X','X','X','X','X','X','X','X','X',
		'X','X','X','X','X','X','X','X','X','X','X','X','X','X','X','X',
		'X','X','X','X','X','X','X','X','X','X','X','X','X','X','X','X',
		'X','X','X','X','X','X','X','X','X','X','X','X','X','X','X','X',
		'X','X','X','X','X','X','X','X','X','X','X','X','X','X','X','X',
		'X','X','X','X','X','X','X','X','X','X','X','X','X','X','X','X',
		'X','X','X','X','X','X','X','X','X','X','X','X','X','X','X','X',
		'X','X','X','X','X','X','X','X','X','X','X','X','X','X','X','X',
		'X','X','X','X','X','X','X','X','X','X','X','X','X','X','X','X',
		'X','X','X','X','X','X','X','X','X','X','X','X','X','X','X','X',
		'X','X','X','X','X','X','X','X','X','X','X','X','X','X','X','X',
		'X','X','X','X','X','X','X','X','X','X','X','X','X','X','X','X',
		'X','X','X','X','X','X','X','X','X','X','X','X','X','X','X','X',
		'X','X','X','X','X','X','X','X','X','X','X','X','X','X','X','X',
		'X','X','X','X','X','X','X','X','X','X','X','X','X','X','X','X',
		'X','X','X','X','X','X','X','X','X','X','X','X','X','X','X','X',
		'X','X','X','X','X','X','X','X','X','X','X','X','X','X','X','X',
		'X','X','X','X','X','X','X','X','X','X','X','X','X','X','X','X',
		'X','X','X','X','X','X','X','X','X','X','X','X','X','X','X','X',
		'X','X','X','X','X','X','X','X','X','X','X','X','X','X','X','X',
		'X','X','X','X','X','X','X','X','X','X','X','X','X','X','X','X',
		'X','X','X','X','X','X','X','X','X','X','X','X','X','X','X','X',
		'X','X','X','X','X','X','X','X','X','X','X','X','X','X','X','X',
		'X','X','X','X','X','X','X','X','X','X','X','X','X','X','X','X',
		'X','X','X','X','X','X','X','X','X','X','X','X','X','X','X','X',
		'X','X','X','X','X','X','X','X','X','X','X','X','X','X','X','X',
		'X','X','X','X','X','X','X','X','X','X','X','X','X','X','X','X',
		'X','X','X','X','X','X','X','X','X','X','X','X','X','X','X','X',
		'X','X','X','X','X','X','X','X','X','X','X','X','X','X','X','X',
		'X','X','X','X','X','X','X','X','X','X','X','X','X','X','X','X',
		'X','X','X','X','X','X','X','X','X','X','X','X','X','X','X','X',
		'X','X','X','X','X','X','X','X','X','X','X','X','X','X','X','X',
		'X','X','X','X','X','X','X','X','X','X','X','X','X','X','X','X',
		'X','X','X','X','X','X','X','X','X','X','X','X','X','X','X','X',
		'X','X','X','X','X','X','X','X','X','X','X','X','X','X','X','X',
		'X','X','X','X','X','X','X','X','X','X','X','X','X','X','X','X',
		'X','X','X','X','X','X','X','X','X','X','X','X','X','X','X','X',
		'X','X','X','X','X','X','X','X','X','X','X','X','X','X','X','X',
		'X','X','X','X','X','X','X','X','X','X','X','X','X','X','X','X',
		'X','X','X','X','X','X','X','X','X','X','X','X','X','X','X','X',
		'X','X','X','X','X','X','X','X','X','X','X','X','X','X','X','X',
		'X','X','X','X','X','X','X','X','X','X','X','X','X','X','X','X',
		'X','X','X','X','X','X','X','X','X','X','X','X','X','X','X','X',
		'X','X','X','X','X','X','X','X','X','X','X','X','X','X','X','X',
		'X','X','X','X','X','X','X','X','X','X','X','X','X','X','X','X',
		'X','X','X','X','X','X','X','X','X','X','X','X','X','X','X','X',
		'X','X','X','X','X','X','X','X','X','X','X','X','X','X','X','X',
		'X','X','X','X','X','X','X','X','X','X','X','X','X','X','X','X',
		'X','X','X','X','X','X','X','X','X','X','X','X','X','X','X','X',
		'X','X','X','X','X','X','X','X','X','X','X','X','X','X','X','X',
		'X','X','X','X','X','X','X','X','X','X','X','X','X','X','X','X',
		'X','X','X','X','X','X','X','X','X','X','X','X','X','X','X','X',
		'X','X','X','X','X','X','X','X','X','X','X','X','X','X','X','X',
		'X','X','X','X','X','X','X','X','X','X','X','X','X','X','X','X',
		'X','X','X','X','X','X','X','X','X','X','X','X','X','X','X','X',
		'X','X','X','X','X','X','X','X','X','X','X','X','X','X','X','X',
		'X','X','X','X','X','X','X','X','X','X','X','X','X','X','X','X',
		'X','X','X','X','X','X','X','X','X','X','X','X','X','X','X','X',
		'X','X','X','X','X','X','X','X','X','X','X','X','X','X','X','X',
		'X','X','X','X','X','X','X','X','X','X','X','X','X','X','X','X',
		'X','X','X','X','X','X','X','X','X','X','X','X','X','X','X','X',
		'X','X','X','X','X','X','X','X','X','X','X','X','X','X','X','X',
		'X','X','X','X','X','X','X','X','X','X','X','X','X','X','X','X',
		'X','X','X','X','X','X','X','X','X','X','X','X','X','X','X','X',
		'X','X','X','X','X','X','X','X','X','X','X','X','X','X','X','X',
		'X','X','X','X','X','X','X','X','X','X','X','X','X','X','X','X',
		'X','X','X','X','X','X','X','X','X','X','X','X','X','X','X','X',
		'X','X','X','X','X','X','X','X','X','X','X','X','X','X','X','X',
		'X','X','X','X','X','X','X','X','X','X','X','X','X','X','X','X',
		'X','X','X','X','X','X','X','X','X','X','X','X','X','X','X','X',
		'X','X','X','X','X','X','X','X','X','X','X','X','X','X','X','X',
		'X','X','X','X','X','X','X','X','X','X','X','X','X','X','X','X',
		'X','X','X','X','X','X','X','X','X','X','X','X','X','X','X','X',
		'X','X','X','X','X','X','X','X','X','X','X','X','X','X','X','X',
		'X','X','X','X','X','X','X','X','X','X','X','X','X','X','X','X',
		'X','X','X','X','X','X','X','X','X','X','X','X','X','X','X','X',
		'X','X','X','X','X','X','X','X','X','X','X','X','X','X','X','X',
		'X','X','X','X','X','X','X','X','X','X','X','X','X','X','X','X'
};

#endif // #ifndef __TRANSLATE_HH
