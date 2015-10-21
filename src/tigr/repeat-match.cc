/* Programmer:  A. Delcher
*     Revised:  12 December 2002
*        File:  MUMmer/repeat-match.cc
*
*  This program identifies maximal exact repeat regions (longer than
*  MIN_MATCH_LEN threshold) in input genome.  Uses both forward and
*  reverse strand.
*/


#include <mummer/tigrinc.hh>


#define  SHOW_TREE 0
#define  USE_EXTRA_FIELDS  0
#define  DEBUG  0


int  Global_Trace = 0;
int  Global_Skip_Ct = 0;
int  Global_Non_Skip_Ct = 0;


#define  NIL  0         // Remove if convert to pointers

const int  DEFAULT_MIN_MATCH_LEN = 20;
const char  DOLLAR_CHAR = '$';
// const char  DONT_KNOW_CHAR = 'N';
const int  MAX_NAME_LEN = 500;
const char  START_CHAR = '%';



typedef  struct node
  {
   int  Lo;
   unsigned  Child : 31;
   unsigned  Child_Is_Leaf : 1;
   int  Sibling : 31;
   unsigned  Sibling_Is_Leaf : 1;
   unsigned  Link;
#if  USE_EXTRA_FIELDS
   int  Depth, ID;
   unsigned  Parent;
#endif
   int  Len : 31;
   unsigned  Should_Skip : 1;
   int  Subtree_Size;
  }  Node_Type, * Node_Ptr;

typedef  struct leaf
  {
   int  Lo;
   int  Sibling : 31;
   unsigned  Sibling_Is_Leaf : 1;
#if  USE_EXTRA_FIELDS
   int  Depth, ID;
   unsigned  Parent;
#endif
   int  Len : 31;
   unsigned  Is_Duplicate : 1;
  }  Leaf_Type, * Leaf_Ptr;


char  * Data;
Leaf_Ptr  Leaf_Array;
Node_Ptr  Node_Array;
int  * Next_Leaf;

int  Curr_ID;
int  Curr_String_ID;
int  Data_Len = 2;
static bool  Exhaustive_Matches = false;
  // Set by -E option; if true then matches are found by exhaustive search
  // For testing purposes
static bool  Forward_Only = false;
  // Set by -f option; if true then matches to reverse complement string
  // are not considered
long int  Genome_Len;
static char  * Input_File_Name = NULL;
  // Name of input fasta file on command line
int  Longest_String = 0;
int  Max_Depth = 0;
int  Min_Match_Len = DEFAULT_MIN_MATCH_LEN;
int  Next_Avail_Node = 1;
int  Num_Strings = 2;
int  String_Separator;
static bool  Tandem_Only = false;
  // Set by -t option to output only tandem repeats
int  Verbose = 0;
  // Set by -V option to do extra tests and/or print debugging output


int  Add_Duplicates
    (int, int, int, int);
int  Add_String
    (int, int);
int  Build_Suffix_Tree
    (int);
static int  New_Find_Child
    (int, char, int &, int &, int &, int);
static int  New_Jump_Down
    (int, int, int, int, int &, int &);
void  List_Matches
    (int, int, int);
int  List_Tree
    (int root, int is_leaf, int parent, int parent_depth);
void  List_Maximal_Matches
    (int, int, int, int);
int  Longest_Prefix_Match
    (char * p, char * q);
void  Mark_Skipable_Nodes
    (int, int, int, int);
int  New_Node
    (void);
static void  Parse_Command_Line
    (int argc, char * argv []);
void  Set_Subtree_Size
    (int, int, int, int);
static int  New_Step_Down
    (int, int, int, int, int &, int, int &, int &);
static void  Usage
    (char * command);
void  Verify_Match
    (int a, int b, int n, int reverse);



int main  (int argc, char * argv [])
 
  {
   FILE  * fp;
   long int  Input_Size;
   int  Root;
   char  Name [MAX_NAME_LEN];

   Parse_Command_Line (argc, argv);

   if  (Verbose > 0)
       {
        fprintf (stderr, "Verbose = %d\n", Verbose);
        printf ("Node size = %d\n", (int) sizeof (Node_Type));
        printf ("Leaf size = %d\n", (int) sizeof (Leaf_Type));
        printf ("sizeof (int) = %d\n", (int) sizeof (int));
        printf ("sizeof (size_t) = %d\n", (int) sizeof (size_t));
       }

   fp = File_Open (Input_File_Name, "r");

   Data = (char *) Safe_malloc (INIT_SIZE);
   Input_Size = INIT_SIZE;

   Read_String (fp, Data, Input_Size, Name, false);
   fclose (fp);

   Genome_Len = strlen (Data + 1);
   Data [0] = START_CHAR;
   Data = (char *) Safe_realloc (Data, 2 * Genome_Len + 4);
   strcpy (Data + Genome_Len + 2, Data + 1);
   Reverse_Complement (Data, Genome_Len + 2, 2 * Genome_Len + 1);
   Data [1 + Genome_Len] = DOLLAR_CHAR;
   Data [2 + 2 * Genome_Len] = DOLLAR_CHAR;
   Data_Len = 3 + 2 * Genome_Len;
   Data [Data_Len] = '\0';

   String_Separator = 1 + Genome_Len;

   Leaf_Array = (Leaf_Ptr) Safe_calloc (Data_Len, sizeof (Leaf_Type));

   Node_Array = (Node_Ptr) Safe_calloc (Data_Len, sizeof (Node_Type));
   Next_Leaf = (int *) Safe_calloc (Data_Len, sizeof (int));

   Curr_String_ID = 0;
   Root = Build_Suffix_Tree (1);

   if  (! Forward_Only)
       {
        Curr_String_ID = 1;
        if  (Add_String (2 + Genome_Len, Root) != 0)
            return  -1;
       }

#if  SHOW_TREE
   printf ("\n\nNodes in Entire Suffix Tree:\n\n");
   List_Tree (Root, false, 0, 0);
#endif

   if(Verbose > 0)
     fprintf (stderr, "Genome Length = %ld   Used %d internal nodes\n",
              Genome_Len, Next_Avail_Node);

   Set_Subtree_Size (Root, false, NIL, 0);

   Mark_Skipable_Nodes (Root, false, NIL, 0);

   if  (Exhaustive_Matches)
       {  // Exhaustive search for matches
        int  i, j, match;

        Data [Genome_Len + 1] = '#';

        printf ("Exhaustive Exact Matches:\n");
        printf ("%9s %10s  %8s\n", "Start1", "Start2", "Length");

        for  (i = 1;  i <= Genome_Len - Min_Match_Len;  i ++)
          for  (j = i + 1;  j <= 1 + Genome_Len - Min_Match_Len;  j ++)
            {
             if  (Data [i - 1] == Data [j - 1])
                 continue;
             match = Longest_Prefix_Match (Data + i, Data + j);
             if  (match >= Min_Match_Len)
                 printf ("%9d %10d%c %8d\n", i, j, ' ', match);
            }

        for  (i = 1;  i <= Genome_Len + 1 - Min_Match_Len;  i ++)
          for  (j = Genome_Len + 2;  j <= 2 * Genome_Len - i;  j ++)
            {
             if  (Data [i - 1] == Data [j - 1])
                 continue;
             match = Longest_Prefix_Match (Data + i, Data + j);
             if  (match >= Min_Match_Len)
                 printf ("%9d %10d%c %8d\n", i, int (2 * Genome_Len + 2 - j), 'r', match);
            }

        Data [Genome_Len + 1] = DOLLAR_CHAR;
       }
     else
       {
        printf ("Long Exact Matches:\n");
        printf ("%9s %10s  %8s\n", "Start1", "Start2", "Length");

        List_Maximal_Matches (Root, false, NIL, 0);
       }

   if  (Verbose > 1)
       {
        fprintf (stderr, "Global_Skip_Ct = %d\n", Global_Skip_Ct);
        fprintf (stderr, "Global_Non_Skip_Ct = %d\n", Global_Non_Skip_Ct);
       }

   return  0;
  }



int  Add_Duplicates  (int Start, int End, int Leaf, int Leaf_Depth)

/* Mark all duplicate occurrences of suffixes in  Data [Start .. End]
*  in suffix tree where the first duplicate is at  Leaf  whose
*  depth is  Leaf_Depth . */

  {
   int  i, j;

   j = Leaf;
   for  (i = Start;  i <= End;  i ++, j++)
     Leaf_Array [j] . Is_Duplicate = true;

   return  0;
  }



int  Add_String  (int Start, int Root)

/* Add all suffixes of string  Data [Start ...]  ending at the
*  first DOLLAR_CHAR encountered to the suffix tree rooted at
*  Root . */

  {
   int  Leaf, End, Last_Parent, Grandparent;
   int  Segment_Len, Segment_Start, Leaf_Len;
   int  Last_Parent_Depth, Link_Depth, Matched, Offset;
   int  New_Place, New_Depth, Made_New_Node, New_Place_Is_Leaf;

   for  (End = Start;  Data [End] != DOLLAR_CHAR;  End ++)
     ;

   Curr_ID = Start;
   Segment_Start = Start;
   Segment_Len = 1 + End - Start;
   New_Place = New_Step_Down (Root, 0, Segment_Start, Segment_Len,
                          Matched, true, New_Place_Is_Leaf, Grandparent);

   if  (Segment_Len == Matched)
       {
        printf ("*** Genome is exact palindrome ***\n");
        return  -1;
       }

   if  (1 + Matched == Segment_Len)
       {
        Offset = 0;
        if  (New_Place_Is_Leaf)
            fprintf (stderr, "ERROR:  Unexpected leaf\n");
          else
            {
             while  (! Node_Array [New_Place] . Child_Is_Leaf)
               {
                New_Place = Node_Array [New_Place] . Child;
                Offset += abs (Node_Array [New_Place] . Len);
               }
             New_Place = Node_Array [New_Place] . Child;
             Offset += abs (Leaf_Array [New_Place] . Len) - 1;
            }
        printf ("String %d is a substring of previous string.\n",
                     Curr_String_ID);
        return  0;
       }

   Leaf = Start;
   Leaf_Array [Leaf] . Lo = Segment_Start + Matched;
   Leaf_Array [Leaf] . Sibling = Node_Array [New_Place] . Child;
   Leaf_Array [Leaf] . Sibling_Is_Leaf
          = Node_Array [New_Place] . Child_Is_Leaf;
   Node_Array [New_Place] . Child = Leaf;
   Node_Array [New_Place] . Child_Is_Leaf = true;
   Leaf_Array [Leaf] . Len = Segment_Len - Matched;
#if  USE_EXTRA_FIELDS
   Leaf_Array [Leaf] . Depth = 1 + End - Start;
   Leaf_Array [Leaf] . ID = Start;
   Leaf_Array [Leaf] . Parent = New_Place;
#endif

   Last_Parent = New_Place;
   Last_Parent_Depth = Matched;

   while  (++ Start <= End)
     {
      Curr_ID ++;
      if  (Node_Array [Last_Parent] . Link != NIL)
          {
           if  (Last_Parent == Root)
               {
                Segment_Start = Start;
                Segment_Len = 1 + End - Start;
                Link_Depth = Last_Parent_Depth;
               }
             else
               {
                Segment_Start = Start + Last_Parent_Depth - 1;
                Segment_Len = 2 + End - Start - Last_Parent_Depth;
                Link_Depth = Last_Parent_Depth - 1;
               }
           New_Place = New_Step_Down (Node_Array [Last_Parent] . Link,
                                  Link_Depth,
                                  Segment_Start, Segment_Len,
                                  Matched, false,
                                  New_Place_Is_Leaf, Grandparent);
           if  (Matched == Segment_Len)
               {
                New_Depth = Link_Depth + Matched;
                return  Add_Duplicates (Start, End, New_Place, New_Depth);
               }

           Leaf = Start;
           Leaf_Array [Leaf] . Lo = Segment_Start + Matched;
           Leaf_Array [Leaf] . Sibling = Node_Array [New_Place] . Child;
           Leaf_Array [Leaf] . Sibling_Is_Leaf
                  = Node_Array [New_Place] . Child_Is_Leaf;
           Node_Array [New_Place] . Child = Leaf;
           Node_Array [New_Place] . Child_Is_Leaf = true;
           Leaf_Array [Leaf] . Len = Segment_Len - Matched;
#if  USE_EXTRA_FIELDS
           Leaf_Array [Leaf] . Depth = 1 + End - Start;
           Leaf_Array [Leaf] . ID = Start;
           Leaf_Array [Leaf] . Parent = New_Place;
#endif

           Last_Parent = New_Place;
           Last_Parent_Depth = Link_Depth + Matched;
          }
        else
          {
           Leaf_Len = Leaf_Array [Leaf] . Len;
           if  (Grandparent == Root)
               {
                Segment_Start = Start;
                Segment_Len = Node_Array [Last_Parent] . Len - 1;
                Link_Depth = 0;
               }
             else
               {
                Segment_Start = Start + Last_Parent_Depth - 1
                                    - Node_Array [Last_Parent] . Len;
                Segment_Len = Node_Array [Last_Parent] . Len;
                Link_Depth = Last_Parent_Depth - 1
                                    - Node_Array [Last_Parent] . Len;
               }
           New_Place = New_Jump_Down (Node_Array [Grandparent] . Link,
                                  Link_Depth,
                                  Segment_Start, Segment_Len,
                                  Made_New_Node, Grandparent);

           if  (Made_New_Node)
               {
                Leaf = Start;
                Leaf_Array [Leaf] . Lo = Segment_Start + Segment_Len;
                Leaf_Array [Leaf] . Sibling = Node_Array [New_Place] . Child;
                Leaf_Array [Leaf] . Sibling_Is_Leaf
                       = Node_Array [New_Place] . Child_Is_Leaf;
                Node_Array [New_Place] . Child = Leaf;
                Node_Array [New_Place] . Child_Is_Leaf = true;
                Leaf_Array [Leaf] . Len = Leaf_Len;
#if  USE_EXTRA_FIELDS
                Leaf_Array [Leaf] . Depth = 1 + End - Start;
                Leaf_Array [Leaf] . ID = Start;
                Leaf_Array [Leaf] . Parent = New_Place;
#endif
                Node_Array [Last_Parent] . Link = New_Place;
                Last_Parent = New_Place;
                Last_Parent_Depth = Link_Depth + Segment_Len;
               }
             else
               {
                Node_Array [Last_Parent] . Link = New_Place;
                Segment_Start += Segment_Len;
                Link_Depth += Segment_Len;
                Segment_Len = Leaf_Len;
                New_Place = New_Step_Down (New_Place, Link_Depth,
                                       Segment_Start, Segment_Len,
                                       Matched, false,
                                       New_Place_Is_Leaf, Grandparent);
                if  (Matched >= Segment_Len)
                    {
                     New_Depth = Link_Depth + Matched;
                     return  Add_Duplicates (Start, End, New_Place, New_Depth);
                    }

                Leaf = Start;
                Leaf_Array [Leaf] . Lo = Segment_Start + Matched;
                Leaf_Array [Leaf] . Sibling
                        = Node_Array [New_Place] . Child;
                Leaf_Array [Leaf] . Sibling_Is_Leaf
                       = Node_Array [New_Place] . Child_Is_Leaf;
                Node_Array [New_Place] . Child = Leaf;
                Node_Array [New_Place] . Child_Is_Leaf = true;
                Leaf_Array [Leaf] . Len = Segment_Len - Matched;
#if  USE_EXTRA_FIELDS
                Leaf_Array [Leaf] . Depth = 1 + End - Start;
                Leaf_Array [Leaf] . ID = Start;
                Leaf_Array [Leaf] . Parent = New_Place;
#endif

                Last_Parent = New_Place;
                Last_Parent_Depth = Link_Depth + Matched;
               }
          }
     }

   return  0;
  }



int  Build_Suffix_Tree  (int Start)

/* Build a suffix tree for the string  Data [Start ...]  ending at
*  the first DOLLAR_CHAR encountered.  Return the subscript of the
*  root of the resulting tree. */

  {
   int  Root, Leaf, End, Last_Parent, Grandparent;
   int  Segment_Len, Segment_Start, Leaf_Len;
   int  Last_Parent_Depth, Link_Depth, Matched;
   int  New_Place, Made_New_Node, New_Place_Is_Leaf;

   for  (End = Start;  Data [End] != DOLLAR_CHAR;  End ++)
     ;

   Curr_ID = Start;
   Root = New_Node ();
   Leaf = Start;

   Node_Array [Root] . Lo = 0;
   Node_Array [Root] . Child = Leaf;
   Node_Array [Root] . Sibling = NIL;
   Node_Array [Root] . Link = Root;
   Node_Array [Root] . Len = 0;
   Node_Array [Root] . Sibling_Is_Leaf = false;
   Node_Array [Root] . Child_Is_Leaf = true;
#if  USE_EXTRA_FIELDS
   Node_Array [Root] . Depth = 0;
   Node_Array [Root] . ID = Curr_ID;
   Node_Array [Root] . Parent = NIL;
#endif

   Leaf_Array [Leaf] . Lo = Start;
   Leaf_Array [Leaf] . Sibling = NIL;
   Leaf_Array [Leaf] . Len = 1 + End - Start;
   Leaf_Array [Leaf] . Sibling_Is_Leaf = false;
#if  USE_EXTRA_FIELDS
   Leaf_Array [Leaf] . Depth = 1 + End - Start;
   Leaf_Array [Leaf] . ID = Start;
   Leaf_Array [Leaf] . Parent = Root;
#endif

   Last_Parent = Root;
   Last_Parent_Depth = 0;

   while  (++ Start <= End)
     {
      if  (Verbose > 1)
          printf ("Build_Suffix_Tree:  Start = %d\n", Start);

      Curr_ID ++;
      if  (Node_Array [Last_Parent] . Link != NIL)
          {
           if  (Last_Parent == Root)
               {
                Segment_Start = Start;
                Segment_Len = 1 + End - Start;
                Link_Depth = Last_Parent_Depth;
               }
             else
               {
                Segment_Start = Start + Last_Parent_Depth - 1;
                Segment_Len = 2 + End - Start - Last_Parent_Depth;
                Link_Depth = Last_Parent_Depth - 1;
               }
           New_Place = New_Step_Down (Node_Array [Last_Parent] . Link,
                                  Link_Depth,
                                  Segment_Start, Segment_Len,
                                  Matched, false,
                                  New_Place_Is_Leaf, Grandparent);
           if  (Matched >= Segment_Len)
               {
                printf ("Ooops:  Suffix can't appear twice.\n");
                exit  (EXIT_FAILURE);
               }
             else
               {
                Leaf = Start;
                Leaf_Array [Leaf] . Lo = Segment_Start + Matched;
                Leaf_Array [Leaf] . Sibling = Node_Array [New_Place] . Child;
                Leaf_Array [Leaf] . Sibling_Is_Leaf
                       = Node_Array [New_Place] . Child_Is_Leaf;
                Node_Array [New_Place] . Child = Leaf;
                Node_Array [New_Place] . Child_Is_Leaf = true;
                Leaf_Array [Leaf] . Len = Segment_Len - Matched;
#if  USE_EXTRA_FIELDS
                Leaf_Array [Leaf] . Depth = 1 + End - Start;
                Leaf_Array [Leaf] . ID = Start;
                Leaf_Array [Leaf] . Parent = New_Place;
#endif
               }
           Last_Parent = New_Place;
           Last_Parent_Depth = Link_Depth + Matched;
          }
        else
          {
           Leaf_Len = Leaf_Array [Leaf] . Len;
           if  (Grandparent == Root)
               {
                Segment_Start = Start;
                Segment_Len = Node_Array [Last_Parent] . Len - 1;
                Link_Depth = 0;
               }
             else
               {
                Segment_Start = Start + Last_Parent_Depth - 1
                                    - Node_Array [Last_Parent] . Len;
                Segment_Len = Node_Array [Last_Parent] . Len;
                Link_Depth = Last_Parent_Depth - 1
                                    - Node_Array [Last_Parent] . Len;
               }

           New_Place = New_Jump_Down (Node_Array [Grandparent] . Link,
                                  Link_Depth,
                                  Segment_Start, Segment_Len,
                                  Made_New_Node, Grandparent);
           if  (Made_New_Node)
               {
                Leaf = Start;
                Leaf_Array [Leaf] . Lo = Segment_Start + Segment_Len;
                Leaf_Array [Leaf] . Sibling = Node_Array [New_Place] . Child;
                Leaf_Array [Leaf] . Sibling_Is_Leaf
                       = Node_Array [New_Place] . Child_Is_Leaf;
                Node_Array [New_Place] . Child = Leaf;
                Node_Array [New_Place] . Child_Is_Leaf = true;
                Leaf_Array [Leaf] . Len = Leaf_Len;
#if  USE_EXTRA_FIELDS
                Leaf_Array [Leaf] . Depth = 1 + End - Start;
                Leaf_Array [Leaf] . ID = Start;
                Leaf_Array [Leaf] . Parent = New_Place;
#endif
                Node_Array [Last_Parent] . Link = New_Place;
                Last_Parent = New_Place;
                Last_Parent_Depth = Link_Depth + Segment_Len;
               }
             else
               {
                Node_Array [Last_Parent] . Link = New_Place;
                Segment_Start += Segment_Len;
                Link_Depth += Segment_Len;
                Segment_Len = Leaf_Len;
                New_Place = New_Step_Down (New_Place, Link_Depth,
                                       Segment_Start, Segment_Len,
                                       Matched, false,
                                       New_Place_Is_Leaf, Grandparent);
                if  (Matched >= Segment_Len)
                    {
                     printf ("Ooops:  Suffix can't appear twice.\n");
                     exit  (EXIT_FAILURE);
                    }
                  else
                    {
                     Leaf = Start;
                     Leaf_Array [Leaf] . Lo = Segment_Start + Matched;
                     Leaf_Array [Leaf] . Sibling
                             = Node_Array [New_Place] . Child;
                     Leaf_Array [Leaf] . Sibling_Is_Leaf
                            = Node_Array [New_Place] . Child_Is_Leaf;
                     Node_Array [New_Place] . Child = Leaf;
                     Node_Array [New_Place] . Child_Is_Leaf = true;
                     Leaf_Array [Leaf] . Len = Segment_Len - Matched;
#if  USE_EXTRA_FIELDS
                     Leaf_Array [Leaf] . Depth = 1 + End - Start;
                     Leaf_Array [Leaf] . ID = Start;
                     Leaf_Array [Leaf] . Parent = New_Place;
#endif
                    }
                Last_Parent = New_Place;
                Last_Parent_Depth = Link_Depth + Matched;
               }
          }
     }

   return  Root;
  }



static int  New_Find_Child
    (int Node, char Ch, int & Is_Leaf, int & Pred, int & Pred_Is_Leaf,
     int Node_Depth)

/* Return the subscript of child of  Node  whose string starts
*  with  Ch  and set  Is_Leaf  to indicate what kind of node
*  that child is.  Set  Pred  to the subscript of the prior
*  sibling of the chosen child and  Pred_Is_Leaf  to indicate
*  whether  Pred  is a leaf or not.   Node_Depth  is the depth of  Node
*  in the tree.
*/

  {
   int  Leaf_Lo, Leaf_Len;
   int  i;
   char  Start_Ch;

   Pred = NIL;
   Pred_Is_Leaf = false;
   i = Node_Array [Node] . Child;
   Is_Leaf = Node_Array [Node] . Child_Is_Leaf;
   if  (Is_Leaf)
       {
        Leaf_Lo = Leaf_Array [i] . Lo;
        Leaf_Len = Leaf_Array [i] . Len;
        if  (Leaf_Len > 0)
            Start_Ch = Data [Leaf_Lo];
          else
            Start_Ch = Complement (Data [Leaf_Lo]);
       }
     else
       {
        if  (Node_Array [i] . Len > 0)
            Start_Ch = Data [Node_Array [i] . Lo];
          else
            Start_Ch = Complement (Data [Node_Array [i] . Lo]);
       }

   while  (i != NIL && Start_Ch != Ch)
     {
      Pred = i;
      Pred_Is_Leaf = Is_Leaf;
      if  (Is_Leaf)
          {
           Is_Leaf = Leaf_Array [i] . Sibling_Is_Leaf;
           i = Leaf_Array [i] . Sibling;
          }
        else
          {
           Is_Leaf = Node_Array [i] . Sibling_Is_Leaf;
           i = Node_Array [i] . Sibling;
          }
      if  (Is_Leaf)
          {
           Leaf_Lo = Leaf_Array [i] . Lo;
           Leaf_Len = Leaf_Array [i] . Len;
           if  (Leaf_Len > 0)
               Start_Ch = Data [Leaf_Lo];
             else
               Start_Ch = Complement (Data [Leaf_Lo]);
          }
        else
          {
           if  (Node_Array [i] . Len > 0)
               Start_Ch = Data [Node_Array [i] . Lo];
             else
               Start_Ch = Complement (Data [Node_Array [i] . Lo]);
          }
     }

   return  i;
  }



static int  New_Jump_Down
    (int Node, int Node_Depth, int Lo, int Len, int & Made_New_Node,
     int & Par_New_Node)

/* Return the subscript of the node that represents the
*  descendant of  Node  that matches  Data [Lo .. Lo + Len - 1] .
*  Check only the first character of each child substring--the
*  rest are assumed to match automatically.
*  Node_Depth  is the depth of  Node  in the suffix tree.
*  If necessary, allocate a new node and return its subscript.
*  Set  Made_New_Node  to TRUE iff a new node was allocated and
*  set  Par_New_Node  to the parent of that new node.
*/

  {
   int  P, Q, P_Is_Leaf, D, Pred, Pred_Is_Leaf;
   int  P_Sib, P_Sib_Is_Leaf;

   Made_New_Node = false;
   if  (Len == 0)
       return  Node;

   if  (Verbose > 1)
       printf ("Jump_Down:  Node = %d  Depth = %d  Lo = %d  Len = %d\n",
               Node, Node_Depth, Lo, Len);

   P = New_Find_Child (Node, Data [Lo], P_Is_Leaf, Pred, Pred_Is_Leaf,
                   Node_Depth);

   while  (P != NIL)
     {
      if  (P_Is_Leaf)
          {
           D = Leaf_Array [P] . Len;
           P_Sib = Leaf_Array [P] . Sibling;
           P_Sib_Is_Leaf = Leaf_Array [P] . Sibling_Is_Leaf;
          }
        else
          {
           D = Node_Array [P] . Len;
           P_Sib = Node_Array [P] . Sibling;
           P_Sib_Is_Leaf = Node_Array [P] . Sibling_Is_Leaf;
          }

      if  (Len == D)
          return  P;

      if  (D < Len)
          {
           Lo += D;
           Len -= D;
           Node = P;
           Node_Depth += D;
           P = New_Find_Child (Node, Data [Lo], P_Is_Leaf,
                           Pred, Pred_Is_Leaf, Node_Depth);
          }
        else
          {
           Q = New_Node ();
           Node_Array [Q] . Lo = Lo;
           Node_Array [Q] . Len = Len;
           Node_Array [Q] . Child = P;
           Node_Array [Q] . Sibling = P_Sib;
           Node_Array [Q] . Sibling_Is_Leaf = P_Sib_Is_Leaf;
           Par_New_Node = Node;
           Node_Array [Q] . Link = NIL;
           Node_Array [Q] . Child_Is_Leaf = P_Is_Leaf;
#if  USE_EXTRA_FIELDS
           Node_Array [Q] . Parent = Node;
           Node_Array [Q] . Depth = Node_Array [Node] . Depth + Len;
           Node_Array [Q] . ID = Curr_ID;
#endif
           if  (Pred == NIL)
               {
                Node_Array [Node] . Child = Q;
                Node_Array [Node] . Child_Is_Leaf = false;
               }
           else if  (Pred_Is_Leaf)
               {
                Leaf_Array [Pred] . Sibling = Q;
                Leaf_Array [Pred] . Sibling_Is_Leaf = false;
               }
             else
               {
                Node_Array [Pred] . Sibling = Q;
                Node_Array [Pred] . Sibling_Is_Leaf = false;
               }

           if  (! P_Is_Leaf)
               {
                Node_Array [P] . Lo += Len;
                Node_Array [P] . Len -= Len;
#if  USE_EXTRA_FIELDS
                Node_Array [P] . Parent = Q;
#endif
                Node_Array [P] . Sibling = NIL;
                Node_Array [P] . Sibling_Is_Leaf = false;
               }
             else
               {
                Leaf_Array [P] . Lo += Len;
                Leaf_Array [P] . Len -= Len;
#if  USE_EXTRA_FIELDS
                Leaf_Array [P] . Parent = Q;
#endif
                Leaf_Array [P] . Sibling = NIL;
                Leaf_Array [P] . Sibling_Is_Leaf = false;
               }
           
           Made_New_Node = true;
           return  Q;
          }

     }

   printf ("Ooops:  Couldn't find appropriate child node\n");
   exit  (EXIT_FAILURE);

   return  0;
  }



void  List_Matches
    (int A, int B, int n)

/* List all maximal matches of length  n   between entries in
*  the list starting at
*  subscript  A  in array  Next_Leaf  and the list starting
*  at subscript  B  in the same array.  Don't list a match
*  if it's already been listed in a different order (because
*  of the reverse complement strand. */

  {
   int  i, j, k, L, R, Reversed;
  
   for  (i = A;  i != NIL;  i = Next_Leaf [i])
     {
      for  (j = B;  j != NIL;  j = Next_Leaf [j])
        {
         if  (Data [i - 1] == Data [j - 1]
                  || Data [i + n] == Data [j + n]
                  || (i > String_Separator && j > String_Separator))
             continue;
         Reversed = false;
         if  (j > String_Separator)
             {
              k = Genome_Len - (j - String_Separator) - n + 2;
              if  (k < i)
                  continue;
              L = i;
              R = k + n - 1;
              Reversed = true;
             }
         else if  (i > String_Separator)
             {
              k = Genome_Len - (i - String_Separator) - n + 2;
              if  (k < j)
                  continue;
              L = j;
              R = k + n - 1;
              Reversed = true;
             }
         else if  (i < j)
             {
              L = i;
              R = j;
             }
           else
             {
              L = j;
              R = i;
             }
         if  (Verbose > 1)
             Verify_Match (L, R, n, Reversed);
         if  (Tandem_Only && L + n < R)
             continue;
         printf ("%9d %10d%c %8d\n", L, R, Reversed ? 'r' : ' ', n);
        }
     }
  }



int  List_Tree  (int Root, int Is_Leaf, int Parent, int Parent_Depth)

//  Show contents of suffix tree rooted at  Root  whose parent's
//  string depth in the suffix tree is  Depth .   Is_Leaf
//  indicates whether  Root  is a leaf node or not.
//   Parent  is the subscript of this node's parent.

  {
   int  i, j;

   while  (Root != NIL)
     if  (Is_Leaf)
         {
          if  (Leaf_Array [Root] . Len > 0)
              {
               j = Leaf_Array [Root] . Lo - Parent_Depth;
               printf ("Leaf %d:  ", j);
               for  (i = 0;  i < Leaf_Array [Root] . Len;  i ++)
                 putchar (Data [Leaf_Array [Root] . Lo + i]);
              }
            else
              {
               j = - (Leaf_Array [Root] . Lo + Parent_Depth);
               printf ("Leaf %d:  ", j);
               for  (i = 0;  i > Leaf_Array [Root] . Len;  i --)
                 putchar (Complement (Data [Leaf_Array [Root] . Lo + i]));
              }
          printf ("  Par = %d", Parent);
#if  USE_EXTRA_FIELDS
          printf ("  Depth = %d", Leaf_Array [Root] . Depth);
#endif
          if  (Leaf_Array [Root] . Is_Duplicate)
              printf ("  Duplicate = %d", int (Root + Genome_Len + 1));
          putchar ('\n');
          Is_Leaf = Leaf_Array [Root] . Sibling_Is_Leaf;
          Root = Leaf_Array [Root] . Sibling;
         }

       else
         {
#if  USE_EXTRA_FIELDS
          printf ("Node %d:  ID %d  ", Root, Node_Array [Root] . ID);
#else
          printf ("Node %d:  ", Root);
#endif
          if  (Node_Array [Root] . Len > 0)
              for  (i = 0;  i < Node_Array [Root] . Len;  i ++)
                putchar (Data [Node_Array [Root] . Lo + i]);
            else
              for  (i = 0;  i > Node_Array [Root] . Len;  i --)
                putchar (Complement (Data [Node_Array [Root] . Lo + i]));
          printf ("  Par = %d", Parent);
#if  USE_EXTRA_FIELDS
          printf ("  Depth = %d", Node_Array [Root] . Depth);
          if  (Node_Array [Root] . Link != NIL)
              printf ("   Link to node %d (ID %d)", Node_Array [Root] . Link,
                      Node_Array [Node_Array [Root] . Link] . ID);
#else
          if  (Node_Array [Root] . Link != NIL)
              printf ("   Link to node %d", Node_Array [Root] . Link);
#endif
          putchar ('\n');
          List_Tree (Node_Array [Root] . Child,
                            Node_Array [Root] . Child_Is_Leaf, Root,
                            Parent_Depth + abs (Node_Array [Root] . Len));
          Is_Leaf = Node_Array [Root] . Sibling_Is_Leaf;
          Root = Node_Array [Root] . Sibling;
         }

   return  0;
  }



void  List_Maximal_Matches  (int Root, int Is_Leaf, int Parent, int Parent_Depth)

/* List substring pairs that occur more than once where the match
*  does not extend either to the left or to the right for
*  the subtree rooted at  Root .  Is_Leaf  indicates
*  whether  Root  is a leaf or internal node.   Parent  is the
*  parent of  Root  and  Parent_Depth  is its string-depth in the
*  suffix tree. */

  {
   int  j, k, Depth, List1, Start, End;
   int  k_Is_Leaf;

   if  (Root == NIL)
       return;

   if  (Is_Leaf)
       {
        // Skip any match here, it's other version at the start of
        // the string will be output instead

        List_Maximal_Matches (Leaf_Array [Root] . Sibling,
               Leaf_Array [Root] . Sibling_Is_Leaf,
               Parent, Parent_Depth);

        return;
       }

   Depth = Parent_Depth + Node_Array [Root] . Len;

   List_Maximal_Matches (Node_Array [Root] . Child,
          Node_Array [Root] . Child_Is_Leaf, Root, Depth);

   List_Maximal_Matches (Node_Array [Root] . Sibling,
          Node_Array [Root] . Sibling_Is_Leaf, Parent, Parent_Depth);

   if  (Verbose > 1)
       printf ("List_Maximal_Matches:  Root = %d  Parent = %d  Parent_Depth = %d\n"
               "   Depth = %d  Subtree_Size = %d\n",
               Root, Parent, Parent_Depth, Depth,
               Node_Array [Root] . Subtree_Size);

if  (Depth >= Min_Match_Len)
    {
     if  (Node_Array [Root] . Should_Skip)
         Global_Skip_Ct ++;
       else
         Global_Non_Skip_Ct ++;
    }
   if  (Depth >= Min_Match_Len
          && ! Node_Array [Root] . Should_Skip)
       {
        Is_Leaf = Node_Array [Root] . Child_Is_Leaf;
        for  (j = Node_Array [Root] . Child;  j != NIL; )
          {
           if  (Verbose > 1)
               printf ("  Child = %d %s\n", j, Is_Leaf ? "leaf" : "node");
           if  (Is_Leaf)
               {
                List1 = j;
                k = Leaf_Array [j] . Sibling;
                k_Is_Leaf = Leaf_Array [j] . Sibling_Is_Leaf;
                Is_Leaf = Leaf_Array [j] . Sibling_Is_Leaf;
                j = Leaf_Array [j] . Sibling;
               }
             else
               {
                List1 = Node_Array [j] . Link;
                k = Node_Array [j] . Sibling;
                k_Is_Leaf = Node_Array [j] . Sibling_Is_Leaf;
                Is_Leaf = Node_Array [j] . Sibling_Is_Leaf;
                j = Node_Array [j] . Sibling;
               }
           while  (k != NIL)
             {
              if  (k_Is_Leaf)
                  {
                   List_Matches (List1, k, Depth);
                   k_Is_Leaf = Leaf_Array [k] . Sibling_Is_Leaf;
                   k = Leaf_Array [k] . Sibling;
                  }
                else
                  {
                   List_Matches (List1, Node_Array [k] . Link, Depth);
                   k_Is_Leaf = Node_Array [k] . Sibling_Is_Leaf;
                   k = Node_Array [k] . Sibling;
                  }
             }
#if  0
           if  (Is_Leaf)
               {
                printf ("  Lf%7d\n", j);
                Is_Leaf = Leaf_Array [j] . Sibling_Is_Leaf;
                j = Leaf_Array [j] . Sibling;
               }
             else
               {
                printf ("  Nd%7d\n", j);
                Is_Leaf = Node_Array [j] . Sibling_Is_Leaf;
                j = Node_Array [j] . Sibling;
               }
#endif
          }
       }

   Start = End = 0;
   Is_Leaf = Node_Array [Root] . Child_Is_Leaf;
   for  (j = Node_Array [Root] . Child;  j != NIL; )
     {
      if  (Is_Leaf)
          {
           if  (Start == 0)
               Start = End = j;
             else
               {
                Next_Leaf [End] = j;
                End = j;
               }
           Is_Leaf = Leaf_Array [j] . Sibling_Is_Leaf;
           j = Leaf_Array [j] . Sibling;
          }
        else
          {
           if  (Start == 0)
               Start = End = Node_Array [j] . Link;
             else
               Next_Leaf [End] = Node_Array [j] . Link;
           while  (Next_Leaf [End] != NIL)
             End = Next_Leaf [End];
           Is_Leaf = Node_Array [j] . Sibling_Is_Leaf;
           j = Node_Array [j] . Sibling;
          }
     }
   Node_Array [Root] . Link = Start;
#if  0
   printf  ("  Leaves: ");
   for  (j = Node_Array [Root] . Link;  j != NIL;  j = Next_Leaf [j])
     printf (" %3d", j);
   printf ("\n");
#endif

   return;
  }



int  Longest_Prefix_Match
    (char * p, char * q)

//  Return the length of the longest common prefix of strings  p
//  and  q .  Assumes they will mismatch before running off the
//  end of either string.

  {
   int  i;

   for  (i = 0;  p [i] == q [i];  i ++)
     ;

   return  i;
  }



void  Mark_Skipable_Nodes
    (int Root, int Is_Leaf, int Parent, int Parent_Depth)

//  Set the  Should_Skip  field of the node that  Root  links to
//  if its  Subtree_Size  is the same as  Root 's,
//  and do the same recursively in  Root 's subtree.
//   Is_Leaf  indicates whether  Root  is a leaf or internal node.
//   Parent  is the parent of  Root  and  Parent_Depth  is its
//  string-depth in the suffix tree.

  {
   int  Link, Depth;

   if  (Root == NIL)
       return;

   if  (Is_Leaf)
       {
        Mark_Skipable_Nodes (Leaf_Array [Root] . Sibling,
               Leaf_Array [Root] . Sibling_Is_Leaf,
               Parent, Parent_Depth);

        return;
       }

   Depth = Parent_Depth + Node_Array [Root] . Len;

   Mark_Skipable_Nodes (Node_Array [Root] . Child, Node_Array [Root] . Child_Is_Leaf,
          Root, Depth);

   Mark_Skipable_Nodes (Node_Array [Root] . Sibling, Node_Array [Root] . Sibling_Is_Leaf,
          Parent, Parent_Depth);

   if  (Verbose > 1)
       printf ("Mark_Skipable_Nodes:  Root = %d  Parent = %d  Parent_Depth = %d\n",
               Root, Parent, Parent_Depth);

   Link = Node_Array [Root] . Link;
   if  (Link != NIL
          && Node_Array [Link] . Subtree_Size == Node_Array [Root] . Subtree_Size)
        Node_Array [Link] . Should_Skip = true;

   // Will mark the root of the entire tree as skippable, but that's
   // what we want

   return;
  }



int  New_Node  (void)

/* Allocate and return the subscript of a new node. */

  {
   return  Next_Avail_Node ++;
  }



static void  Parse_Command_Line
    (int argc, char * argv [])

//  Get options and parameters from command line with  argc
//  arguments in  argv [0 .. (argc - 1)] .

  {
   bool  errflg = false;
   char  * p;
   int  ch;

   optarg = NULL;

   while  (! errflg && ((ch = getopt (argc, argv, "Efn:tV:")) != EOF))
     switch  (ch)
       {
        case  'E' :
          Exhaustive_Matches = true;
          break;

        case  'f' :
          Forward_Only = true;
          break;

        case  'n' :
          Min_Match_Len = (int) strtol (optarg, & p, 10);          
          if  (p == optarg || Min_Match_Len < 1)
              {
               fprintf (stderr, "ERROR:  Illegal min match length \"%s\"\n",
                        optarg);
               exit (-1);
              }
          break;

        case  't' :
          Forward_Only = true;
          Tandem_Only = true;
          break;

        case  'V' :
          Verbose = (int) strtol (optarg, & p, 10);          
          if  (p == optarg)
              {
               fprintf (stderr, "ERROR:  Illegal verbose number \"%s\"\n",
                        optarg);
               exit (-1);
              }
          break;

        case  '?' :
          fprintf (stderr, "Unrecognized option -%c\n", optopt);

        default :
          errflg = true;
       }

   if  (errflg || optind != argc - 1)
       {
        Usage (argv [0]);
        exit (EXIT_FAILURE);
       }

   Input_File_Name = argv [optind ++];

   return;
  }



void  Set_Subtree_Size
    (int Root, int Is_Leaf, int Parent, int Parent_Depth)

//  Set the  Subtree_Size  field of  Root  and all its descendants
//  to the number of leaves in its subtree.
//   Is_Leaf  indicates whether  Root  is a leaf or internal node.
//   Parent  is the parent of  Root  and  Parent_Depth  is its
//  string-depth in the suffix tree.

  {
   int Depth;
   //   int Start;

   if  (Root == NIL)
       return;

   if  (Is_Leaf)
       {
         //        Start = Leaf_Array [Root] . Lo - Parent_Depth;
        if  (Leaf_Array [Root] . Is_Duplicate)
            Node_Array [Parent] . Subtree_Size ++;
          else
            Node_Array [Parent] . Subtree_Size += 2;

        Set_Subtree_Size (Leaf_Array [Root] . Sibling,
               Leaf_Array [Root] . Sibling_Is_Leaf,
               Parent, Parent_Depth);

        return;
       }

   Depth = Parent_Depth + Node_Array [Root] . Len;

   Set_Subtree_Size (Node_Array [Root] . Child, Node_Array [Root] . Child_Is_Leaf,
          Root, Depth);

   Set_Subtree_Size (Node_Array [Root] . Sibling, Node_Array [Root] . Sibling_Is_Leaf,
          Parent, Parent_Depth);

   if  (Verbose > 1)
       printf ("Set_Subtree_Size:  Root = %d  Parent = %d  Parent_Depth = %d\n",
               Root, Parent, Parent_Depth);

   if  (Parent != NIL)
       Node_Array [Parent] . Subtree_Size += Node_Array [Root] . Subtree_Size;

   return;
  }



static int  New_Step_Down
    (int Node, int Node_Depth, int Lo, int Len, int & Depth, int Doing_Prefix,
     int & Return_Is_Leaf, int & Par_New_Node)

/* Return the subscript of the node that represents the lowest
*  descendant of  Node  that matches  Data [Lo .. Lo + Len - 1] .
*  If necessary, allocate a new node and return its subscript,
*  unless  Doing_Prefix .  In that case do not allocate a new node
*  if  1 + Depth == Len  but return the child of where the new
*  node would be allocated.  Set  Depth  to the number of characters
*  successfully matched from  Node  down to the returned  Node .
*  Node_Depth  is the depth of  Node  in the suffix tree.
*  Set  Return_Is_Leaf  to indicate the status of the returned
*  subscript.  If a new node is allocated, set  Par_New_Node
*  to its parent.
*/

  {
   int  P, Q, P_Is_Leaf, i, j, D, Pred, Pred_Is_Leaf;
   int  P_Sib, P_Sib_Is_Leaf;

   Depth = 0;

   if  (Len == 0)
       return  NIL;

   P = New_Find_Child (Node, Data [Lo], P_Is_Leaf, Pred, Pred_Is_Leaf, Node_Depth);

   while  (P != NIL)
     {
      if  (P_Is_Leaf)
          {
           i = Leaf_Array [P] . Lo;
           D = Leaf_Array [P] . Len;
           P_Sib = Leaf_Array [P] . Sibling;
           P_Sib_Is_Leaf = Leaf_Array [P] . Sibling_Is_Leaf;
          }
        else
          {
           i = Node_Array [P] . Lo;
           D = Node_Array [P] . Len;
           P_Sib = Node_Array [P] . Sibling;
           P_Sib_Is_Leaf = Node_Array [P] . Sibling_Is_Leaf;
          }
      for  (j = 1;  j < Len && j < D
                      && Data [i + j] == Data [Lo + j];
                 j ++)
        ;
      Depth += j;

      if  (j < D)
          {
           if  (Doing_Prefix && 1 + j == Len)
               {
                Return_Is_Leaf = P_Is_Leaf;
                return  P;
               }
           Q = New_Node ();
           Node_Array [Q] . Lo = Lo;
           Node_Array [Q] . Len = j;
           Node_Array [Q] . Child = P;
           Node_Array [Q] . Sibling = P_Sib;
           Node_Array [Q] . Sibling_Is_Leaf = P_Sib_Is_Leaf;
           Par_New_Node = Node;
           Node_Array [Q] . Link = NIL;
           Node_Array [Q] . Child_Is_Leaf = P_Is_Leaf;
#if  USE_EXTRA_FIELDS
           Node_Array [Q] . Parent = Node;
           Node_Array [Q] . Depth = Node_Array [Node] . Depth + j;
           Node_Array [Q] . ID = Curr_ID;
#endif
           if  (Pred == NIL)
               {
                Node_Array [Node] . Child = Q;
                Node_Array [Node] . Child_Is_Leaf = false;
               }
           else if  (Pred_Is_Leaf)
               {
                Leaf_Array [Pred] . Sibling = Q;
                Leaf_Array [Pred] . Sibling_Is_Leaf = false;
               }
             else
               {
                Node_Array [Pred] . Sibling = Q;
                Node_Array [Pred] . Sibling_Is_Leaf = false;
               }

           if  (! P_Is_Leaf)
               {
                Node_Array [P] . Lo += j;
                Node_Array [P] . Len -= j;
#if  USE_EXTRA_FIELDS
                Node_Array [P] . Parent = Q;
#endif
                Node_Array [P] . Sibling = NIL;
                Node_Array [P] . Sibling_Is_Leaf = false;
               }
             else
               {
                Leaf_Array [P] . Lo += j;
                Leaf_Array [P] . Len -= j;
#if  USE_EXTRA_FIELDS
                Leaf_Array [P] . Parent = Q;
#endif
                Leaf_Array [P] . Sibling = NIL;
                Leaf_Array [P] . Sibling_Is_Leaf = false;
               }
           
           Return_Is_Leaf = false;
           return  Q;
          }

      if  (P_Is_Leaf || j == Len)
          {
//           Return_Is_Leaf = true;
           Return_Is_Leaf = P_Is_Leaf;
           return  P;
          }

      Lo += j;
      Len -= j;
      Node = P;
      Node_Depth += j;
      P = New_Find_Child (Node, Data [Lo], P_Is_Leaf,
                      Pred, Pred_Is_Leaf, Node_Depth);
     }

   Return_Is_Leaf = false;
   return  Node;
  }



static void  Usage
    (char * command)

//  Print to stderr description of options and command line for
//  this program.   command  is the command that was used to
//  invoke it.

  {
   fprintf (stderr,
           "USAGE:  %s  [options]  <genome-file>\n"
           "\n"
           "Find all maximal exact matches in <genome-file>\n"
           "\n"
           "Options:\n"
           " -E    Use exhaustive (slow) search to find matches\n"
           " -f    Forward strand only, don't use reverse complement\n"
           " -n #  Set minimum exact match length to #\n"
           " -t    Only output tandem repeats\n"
           " -V #  Set level of verbose (debugging) printing to #\n"
           "\n",
           command);

   return;
  }



void  Verify_Match
    (int a, int b, int n, int reverse)

//  Verify that the match of length  n  starting at  Data [a]  and
//  Data [b]  is a maximal match.  If  reverse  is true then the
//  string at  Data [b]  goes in the reverse direction.
//  Assume  Data  is padded at the ends to prevent errors
//  caused by falling off the ends.

  {
   int  i;

   if  (reverse)
       {
        if  (Data [a - 1] == Complement (Data [b + 1]))
            printf ("Verify_Match:  a=%d  b=%d  rev=%c  preceding chars match\n",
                    a, b, 't');
        if  (Data [a + n] == Complement (Data [b - n]))
            printf ("Verify_Match:  a=%d  b=%d  rev=%c  following chars match\n",
                    a, b, 't');
        for  (i = 0;  i < n;  i ++)
          if  (Data [a + i] != Complement (Data [b - i]))
              {
               printf ("Verify_Match:  a=%d  b=%d  rev=%c  mismatch at %d[%c]:%d[%c]\n",
                       a, b, 't', a + i, Data [a + i], b - i, Complement (Data [b - i]));
               exit (EXIT_FAILURE);
              }
       }
     else
       {
        if  (Data [a - 1] == Data [b - 1])
            printf ("Verify_Match:  a=%d  b=%d  rev=%c  preceding chars match\n",
                    a, b, 'f');
        if  (Data [a + n] == Data [b + n])
            printf ("Verify_Match:  a=%d  b=%d  rev=%c  following chars match\n",
                    a, b, 'f');
        for  (i = 0;  i < n;  i ++)
          if  (Data [a + i] != Data [b + i])
              {
               printf ("Verify_Match:  a=%d  b=%d  rev=%c  mismatch at %d[%c]:%d[%c]\n",
                       a, b, 'f', a + i, Data [a + i], b + i, Data [b + i]);
               exit (EXIT_FAILURE);
              }
       }

   return;
  }
