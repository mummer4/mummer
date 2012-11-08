# Usage:  awk -f tandem-repeat.awk
#   Outputs tandem repeat regions based on repeat matches found
#   by  repeat-match  program.  That program should be run with
#   the  -t  option (for tandem repeats) or at least the -f
#   option (fo forward strand only), and the output
#   sorted by first and then second column (with the first two
#   header lines removed).

BEGIN   {
         printf "%8s %8s %8s %10s\n", "Start", "Extent", "UnitLen", "Copies";
        }

        {
         if  ($1 + $3 < $2)
             next;
         if  ($1 == prev)
             next;
         start = $1;
         extent = $2 + $3 - $1;
         unitlen = $2 - $1;
         printf "%8d %8d %8d %10.1f\n", start, extent, unitlen, extent / unitlen;
         prev = $1;
        }
