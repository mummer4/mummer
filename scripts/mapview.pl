#!@PERL@

use lib "@LIB_DIR@";
use Foundation;

my $SCRIPT_DIR = "@LIB_DIR@";


my $VERSION_INFO = q~
mapview version 1.01
    ~;


my $HELP_INFO = q~
  USAGE: mapview  [options]  <coords file>  [UTR coords]  [CDS coords]

  DESCRIPTION:
    mapview is a utility program for displaying sequence alignments as
    provided by MUMmer, NUCmer, PROmer or Mgaps. mapview takes the output of
    show-coords and converts it to a FIG, PDF or PS file for visual analysis.
    It can also break the output into multiple files for easier viewing and
    printing.

  MANDATORY:
    coords file    The output of 'show-coords -rl[k]' or 'mgaps'

  OPTIONS:
    UTR coords      UTR coordinate file in GFF format
    CDS coords      CDS coordinate file in GFF format

    -d|maxdist      Set the maximum base-pair distance between linked matches
                    (default 50000)
    -f|format       Set the output format to 'pdf', 'ps' or 'fig'
                    (default 'fig')
    -h
    --help          Display help information and exit
    -m|mag          Set the magnification at which the figure is rendered,
                    this is an option for fig2dev which is used to generate
                    the PDF and PS files (default 1.0)
    -n|num          Set the number of output files used to partition the
                    output, this is to avoid generating files that are too
                    large to display (default 10)
    -p|prefix       Set the output file prefix
                    (default "PROMER_graph or NUCMER_graph")
    -v
    --verbose       Verbose logging of the processed files
    -V
    --version       Display the version information and exit
    -x1 coord       Set the lower coordinate bound of the display
    -x2 coord       Set the upper coordinate bound of the display
    -g|ref          If the input file is provided by 'mgaps', set the
                    reference sequence ID (as it appears in the first column
                    of the UTR/CDS coords file)
    -I              Display the name of query sequences
    -Ir             Display the name of reference genes
    ~;


my $USAGE_INFO = q~
  USAGE:  mapview  [options]  <coords file>  [UTR coords]  [CDS coords]
    ~;


my @DEPEND_INFO =
    (
     "fig2dev",
     "$SCRIPT_DIR/Foundation.pm"
     );

my $err_gff = q~
  ERROR in the input files ! The reference seq ID can't be found in GFF files !
      The first column in the GFF file should be the ID of the reference seq. 
      The alignments file should provide the same info in the column before the last one.
      
      Here are some example records for the GFF file:
     
      gnl|FlyBase|X	Dmel3	initial-exon	2155	2413	.	-	.	X_CG3038.1
      gnl|FlyBase|X	Dmel3	last-exon	1182	2077	.	-	.	X_CG3038.1
      ...
      The fields are :
      <seq_ID> <source> <exon type> <start> <end> <score> <strand> <frame> <gene_name>
    ~;

my $tigr;
my $err;

my $alignm;  
my $futr;
my $fcds;

#-- Initialize TIGR::Foundation
$tigr = new TIGR::Foundation;
if ( !defined ($tigr) ) {
    print (STDERR "ERROR: TIGR::Foundation could not be initialized");
    exit (1);
}
        
#-- Set help and usage information
$tigr->setHelpInfo ($HELP_INFO);
$tigr->setUsageInfo ($USAGE_INFO);
$tigr->setVersionInfo ($VERSION_INFO);
$tigr->addDependInfo (@DEPEND_INFO);

$err = $tigr->TIGR_GetOptions
    (
     "d|maxdist=i" => \$match_dist,
     "f|format=s" => \$format,
     "m|mag=f" => \$magn,
     "n|num=i" => \$noOutfiles,
     "p|prefix=s" => \$outfilename,
     "x1=i" => \$x1win,
     "x2=i" => \$x2win,
     "v|verbose" => \$verb,
     "g|ref=s" => \$Mgaps,
     "I" => \$printIDconting,
     "Ir" => \$printIDgenes
     );

if ( $err == 0  ||  scalar(@ARGV) < 1  || scalar(@ARGV) > 3 ) {
    $tigr->printUsageInfo( );
    print (STDERR "Try '$0 -h' for more information.\n");
    exit (1);
}

($alignm,$futr,$fcds)=@ARGV; 

if ((substr($x1win,0,1) eq '-') || (substr($x2win,0,1) eq '-')){
    print "ERROR2 : coords x1,x2 should be positive integers !!\n";
    $info=1;
}
if ($x1win > $x2win) {
    print "ERROR3 : wrong range coords : x1 >= x2 !!!\n";
    $info=1; 
}
if ($Mgaps){ #formating the mgaps output to be similar with show-coords output
   format_mgaps();
}

if  (!$format){$format="fig";}
if (($x1win) and ($x2win)){  $startfind=0; }
else{  $startfind=1; }
$endfind=0; 
if (!$noOutfiles) {$noOutfiles=10;} 
if (!$match_dist){$match_dist=50000;}

#```````init colors````````````````````````````
$color{"2"}=27;#dark pink 5utr
$color{"3"}=2;#green  ex
$color{"4"}=1;#blue 3utr
#``````````````````````
@linkcolors=(31,14,11);
#``````````````````````````````````````````````
open(F,$alignm);
<F>;
$prog=<F>; chomp($prog);
<F>;
$_=<F>;
@a=(m/\s+(\||\[.+?\])/g) ;
for ($ind=0;$ind<=$#a;$ind++){
   if ($a[$ind] eq "[S1]") {
     $ind_s1=$ind;     
   }
   elsif ($a[$ind] eq "[E1]"){
     $ind_e1=$ind;    
   }
  elsif ($a[$ind] eq "[% IDY]"){
     $ind_pidy=$ind;  
   }
   elsif ($a[$ind] eq "[LEN R]"){
     $ind_lenchr=$ind;   
   }
#   elsif ($a[$ind] eq "[TAGS]"){    
#    $ind_tags=$ind+1; #there are two columns for this header col
#   }   
}
<F>;$mref=-1; 
while (<F>){
   chomp;
   @a=split;                                 
   if (!exists $hRefContigId{$a[-2]}) {       # print  $a[-2]."\n";               
     $hRefContigId{$a[-2]}=$a[$ind_lenchr];
     $lenrefseqs+=$a[$ind_lenchr];  
     $mref++;                               
   }
}
$nobpinfile=int($lenrefseqs/$noOutfiles);         

close(F);

#``````````````````````


if (@ARGV > 1) {
  get_cds_ends();              
  get_utrcds_info();
  test_overlap();
} 
elsif(!$mref) {     
  $fileno=$noOutfiles;         
  $startcoord=0;$endcoord=0;
  for ($i=0;$i<$fileno;$i++) {            
   $endcoord=$startcoord+$nobpinfile-1;
   $endcoord=$lenrefseqs if ($endcoord>$lenrefseqs);
   $file[$i]="$startcoord $endcoord";  
   $startcoord=$endcoord+1;
  }
}

    
$Yorig=3000;
$YdistPID=2000;
$yscale=$YdistPID/50;
$Xscale=14.5;
$gap=800;
#$maxfiles = ($fileno < 10) ? $fileno : 10;
#---------------------------------

if (!$mref){ 
  for($i=0; $i < $fileno; $i++) {
    $nrf=$i;
    set_output_fname();  
    ($startcoord,$endcoord)=split(/\s+/,$file[$i]);  
     open(O,">$procfile".$nrf.".fig");
     print_header(); 
     print $procfile.$nrf.".fig\t range : $startcoord\t$endcoord \n" if ($verb && ($format eq "fig"));
    
     $xs=0;
     $xe=int(($endcoord-$startcoord+1)/$Xscale);
    #$xs=200;
    #$xe=$xs+int(($endcoord-$startcoord+1)/$Xscale);
     print_grid($xs,$xe,$startcoord,$endcoord);
   
     $tmpIdQrycontig="";
     $linkcolor=$linkcolors[0];               
 
     open(F,$alignm);
     <F>;<F>;<F>;<F>;<F>;
     while(<F>) {                             
	chomp;
	@a=split;                            
	if($a[$ind_s1] > $endcoord) { last;}
	if($a[$ind_s1]<$startcoord && $a[$ind_e1] > $startcoord ) { $a[$ind_s1]=$startcoord;}
	if($a[$ind_s1] < $endcoord && $endcoord < $a[$ind_e1]) { $a[$ind_e1]=$endcoord;}
	
	if($a[$ind_s1]>=$startcoord && $a[$ind_e1]<=$endcoord) {
          $x1=int(($a[$ind_s1]-$startcoord)/$Xscale);#
          $x2=int(($a[$ind_e1]-$startcoord)/$Xscale);# 
          print_align($x1,$x2);
	}
     } 
     close(F);
     %hQrycontig=();
     print_genes() if ($futr);
     print_legend();
     close(O);
     change_file_format() if ($format ne "fig");   
  }
}
elsif($mref){#multiple ref seqs
  set_output_fname();
  $tmpIdQrycontig="";
  $linkcolor=$linkcolors[0];
  $startdrawX=0;
  $proclen=0;
  $first=1;
  $nrf=0;
  open(F,$alignm);
  <F>;<F>;<F>;<F>;<F>;
  while(<F>) {
    chomp;
    @a=split;
    if ($a[-2] ne $tmpcontig){
       %hQrycontig=();
       $tmpcontig=$a[-2];
       if ($first){
         $first=0;
         $nrf++;
         open(O,">$procfile".$nrf.".fig");
         print_header();
         print $procfile.$nrf.".fig"."\n" if ($verb && ($format eq "fig"));
         $len=$hRefContigId{$a[-2]};
       }
       else {
         $startdrawX+=int($len/$Xscale)+$gap;
         $len=$hRefContigId{$a[-2]};
         if (($proclen+$len>$nobpinfile) and ($proclen != 0)){
            print_legend();
            close(O);
            change_file_format() if ($format ne "fig"); 
            $nrf++;         
            open(O,">$procfile".$nrf.".fig");
            print_header();   
            print "\n".$procfile.$nrf.".fig"."\n" if ($verb && ($format eq "fig"));
            $proclen=0;
            $startdrawX=0;
         }
       }
       $xs=$startdrawX;
       $xe=$startdrawX+int($len/$Xscale);
       print_grid($xs,$xe,0,$len);
       #print genes from %geneinfo for contig
       print $a[-2]."\t".$hRefContigId{$a[-2]}."\n";
       print_genes_mr() if ($futr);
       $proclen+=$len;      
    }#end if new contig
    $x1=$startdrawX+int($a[$ind_s1]/$Xscale);
    $x2=$startdrawX+int($a[$ind_e1]/$Xscale);   
    print_align($x1,$x2);
  }
  print_legend();
  close(O);
  change_file_format() if ($format ne "fig");

  close(F);
}
#*******************************************************************************
#*******************************************************************************
sub set_output_fname{
  
  if (!$outfilename) {$procfile=$prog."_graph"."_";}
  else {$procfile=$outfilename."_";}
  
  if ($format ne "fig"){ 
     $procfile="tmp".$procfile;
  }
}
#*********************************
sub get_cds_ends{
#3.  print "create \%hcds_ends...\n";
$testGffFormat=0;
 open(F,"<".$fcds);#|| die "can't open \" $fcds cds \" file !";
 while(<F>) {
   chomp;
   if($_) { 
     @a=split;                          
     if (exists $hRefContigId{$a[0]}){#record if at least one of the ref id is the same in GFF and Align files
       $testGffFormat++;
     }
     $genename=$a[8];                                   
     if ($genename ne $tmpname){
        if ($sign eq "+"){ $hcds_ends{$tmpname} = "$cds5 $cds3";}
        elsif ($sign eq "-"){ $hcds_ends{$tmpname} = "$cds3 $cds5";} 
        $tmpname=$genename; 
        $sign=$a[6]; 
     }
     if($sign eq "-") {
	$temp=$a[3]; 
	$a[3]=$a[4];
	$a[4]=$temp;
     }
     if ($a[2] eq "single-exon"){ 
        $cds5=$a[3];
        $cds3=$a[4];
     }
     elsif ($a[2] eq "initial-exon"){ 
        $cds5=$a[3];
     }
     elsif ($a[2] eq "last-exon"){
        $cds3=$a[4];
     }   
   }
 }
 if ($sign eq "+"){ $hcds_ends{$tmpname} = "$cds5 $cds3";}
 elsif ($sign eq "-"){$hcds_ends{$tmpname} = "$cds3 $cds5";} 
                                                           
 test_formatGFF(); 

 #foreach $k ( keys %hcds_ends){
 #  print "cds_ends: ".$k."\t"."\n"; 
 #} exit;
 
}

#*********************************
sub test_formatGFF{
 if ($testGffFormat==0){         
   print (STDERR  "$err_gff \n"); 
   exit (1);
 }

}

#*********************************
sub get_utrcds_info{
#test for gene overlap $geneinfo{gene_name}->[0]=level,stock gene 5'3' utr ends,
#determina %geneinfo{id gene}->5utr,ex,3utr  
#and @file 
# get_gene_ends();
                $testGffFormat=0;
open(F,"<".$futr);# || die "can't open \" $futr utr \" file !";
while(<F>) {
  chomp;
  if($_) { 
     @a=split;
     #get gene ends-utr
     $genename=$a[8];                      
     if (exists $hRefContigId{$a[0]}){ #check if [align file(col before the last one)] = [GFF file(col 1)]
       $testGffFormat++;
     }
           
     if ($genename ne $tmpname){
        if ($sign eq "+"){ $utr_ends{$tmpname} = "$utr5 $utr3";}
        elsif ($sign eq "-"){ $utr_ends{$tmpname} = "$utr3 $utr5";} 
       
        if ($tmpgene) {# for the distinct_utr_cds
          get_utrcds_ends() ;        
        }
        $tmpgene="$a[3] $a[4]";# 
        
        $tmpname=$genename; 
        $sign=$a[6];  
        $hContig_genes{$a[0]}.=" ".$a[8];  # for multiple ref. seqs  
     }
     else{# for the distinct_utr_cds
        if($sign eq "-") {
	$tmpgene="$a[3] $a[4];".$tmpgene;
        }
        else { $tmpgene.=";$a[3] $a[4]"; }
     }#
     
     if($sign eq "-") {
	$temp=$a[3]; 
	$a[3]=$a[4];
	$a[4]=$temp;
     }
     if ($a[2] eq "single-exon"){  
        $utr5=$a[3];
        $utr3=$a[4];
     }
     elsif ($a[2] eq "initial-exon"){ 
        $utr5=$a[3];
     }
     elsif ($a[2] eq "last-exon"){
        $utr3=$a[4];
     }
    
    #init gene info (level) 
    $geneinfo{$a[8]}->[0]=1 if (!exists $geneinfo{$a[8]});                 
        
  }
}
 # for the distinct_utr_cds
 get_utrcds_ends();  
 %cds_ends=();# 

 if ($sign eq "+"){$utr_ends{$tmpname}="$utr5 $utr3"; }
 elsif ($sign eq "-"){$utr_ends{$tmpname}="$utr3 $utr5"; } 
 $hContig_genes{$a[0]}.=" ".$a[8];

 close(F);
 test_formatGFF();
 
}


#****************************************************************
sub get_utrcds_ends{
 $u5="";$ex="";$u3="";
 
 if ($fcds eq $futr) {
   $ex=$tmpgene;
 } 
 else {  
   @ex=split(";",$tmpgene); 
   @cds=split(" ",$hcds_ends{$tmpname});
 
   for($i=0;$i<=$#ex;$i++){
      @coord=split(" ",$ex[$i]);
    
      if ($cds[0]>$coord[0]){
        if ($cds[0]>$coord[1]){
           $u5.="$coord[0] $coord[1];";      
        }
        else{
           $u5.= "$coord[0] "; $u5.=$cds[0]-1 .";" ; #?
           if ($cds[1]<$coord[1]){
             $ex.="$cds[0] $cds[1];";
             $u3.=$cds[1]+1 ." $coord[1];";
           }
           else{ $ex.="$cds[0] $coord[1];";}       
        }    
      }
      else {
        if ($cds[1]>$coord[0]){
          if ($cds[1]>$coord[1]){ 
            $ex.="$coord[0] $coord[1];";       
          }
          else{
            $ex.="$coord[0] $cds[1];";
            $u3.=$cds[1]+1 ." $coord[1];";
          }     
        }
        else { $u3.="$coord[0] $coord[1];";}    
      }    
   }
   chop($u5, $ex, $u3);
 }
 $geneinfo{$tmpname}->[1]=$sign;
 $geneinfo{$tmpname}->[2]=$u5;
 $geneinfo{$tmpname}->[3]=$ex;
 $geneinfo{$tmpname}->[4]=$u3;
}
#*********************************
sub test_overlap{

 if (!$mref){
   $fileno=0;###
   $endcoord=0;###
 }
 foreach $kcontgid (sort keys %hContig_genes){                       
      @allgenes=split(/\s+/,$hContig_genes{$kcontgid});       
      for ($i=1;$i<=$#allgenes;$i++){
         @g1=split (" ", $utr_ends{$allgenes[$i]});    
         $Utr5End{$allgenes[$i]}=$g1[0]; ###
    
         for ($j=$i+1;$j<=$#allgenes;$j++){ #comparing with the rest of the genes
            @g2=split (" ", $utr_ends{$allgenes[$j]});   
            #if the genes are overpaling and they have the same level ,the second gene is liflet to the next level 
            if ( (($g2[0]>=$g1[0]) and ($g2[0]<=$g1[1])) or (($g2[1]>=$g1[0]) and ($g2[1]<=$g1[1])) or (($g1[0]>=$g2[0]) and ($g1[0]<=$g2[1])) ){
              if ($geneinfo{$allgenes[$i]}->[0] == $geneinfo{$allgenes[$j]}->[0]){
                $geneinfo{$allgenes[$j]}->[0]=$geneinfo{$allgenes[$i]}->[0] + 1 ;
              }
            }   
         }                        
         SetTheRangeForEachFile() if ((!$endfind) and (!$mref)); ###         
      }                                              
  }
  $file[$fileno++]="$startcoord $endcoord" if (!$mref);###                        
  %utr_ends=();###  
 
}  

#**************************************
sub SetTheRangeForEachFile{      
   $currstart=$g1[0];
   $currend=$g1[1];
   #---test range ends intersection
   if (!$startfind) {
     if (($x1win <= $currstart) || ($x1win <= $currend)){
         $currstart = $x1win;
         $startfind = 1;
      }  
   }
   if ( $startfind && $x1win && $x2win){
       if (($x2win <= $currstart) || ($x2win <= $currend) ){
           $currend = $x2win;
           $endfind = 1; 
        } 
   }
#--------------------
   if ($startfind) { 
       if(!$endcoord) {
          #$startcoord=0;
          $startcoord = $x1win ? $x1win : 0;
	  $endcoord=$currend;
       }
       else {
         if($currend > $endcoord) {
	    if($currend-$startcoord < $nobpinfile) {
               $endcoord=$currend; 
	    }
	    else {
	       $file[$fileno++]="$startcoord $endcoord";
	       $startcoord=$endcoord+1;
	       $endcoord=$currend;
	    }
	  }
       }
   }#if startfind        
}
#*********************************
sub print_header{
  print O "#FIG 3.2\nLandscape\nCenter\nInches\nLetter  \n100.00\nMultiple\n-2\n1200 2\n";
}
#*********************************
sub print_align{
  my ($x1,$x2)=@_;
 
  $a[$ind_pidy]=50 if ($a[$ind_pidy]<50);
  $a[$ind_pidy]=int($a[$ind_pidy]);
  if ($Mgaps){
     $y=$Yorig+250+$YdistPID-$yscale*2;
     if($a[$#a]=~/rev$/){$y-=25*$yscale;}      
  }
  else{
    $y=$Yorig+250+$YdistPID-$yscale*($a[$ind_pidy]-50);
  }
  if($x1==$x2) { $x2++;}
     #draw the line between matches. is dif color for each contig
  if ($a[$#a] eq $tmpIdQrycontig) {               
     print_connections($hQrycontig{$tmpIdQrycontig}->[1], $x1,$y);
  }
  else{#new contig
   #remember the start coord for printing the id alignments
   if ($printIDconting){                                            
      if  ( $x1 - $XlastPrint > 400 ) { 
         print O "4 0 0 5 0 0 8 0.0000 4 90 270 ";
         printf O ("\t%.0f %.0f ",$x1,$y); 
         print O  $a[$#a], "\\001\n";
         $XlastPrint=$x1; $YlastPrint=$y;
      }
   }
   # 
     ##if it was seen before,but interrupted by another contig   
     if ((exists $hQrycontig{$a[$#a]}) and ($a[$ind_s1]-$hQrycontig{$a[$#a]}->[2] < $match_dist )) { 
        $linkcolor=$hQrycontig{$a[$#a]}->[0];
        print_connections($hQrycontig{$a[$#a]}->[1], $x1,$y);
     }
     else{
         #change the link color 
         unshift(@linkcolors, pop(@linkcolors)); 
         $linkcolor=$linkcolors[0];
         $hQrycontig{$a[$#a]}->[0] = $linkcolor;
     }
   }
   $tmpIdQrycontig=$a[$#a];
   $hQrycontig{$tmpIdQrycontig}->[1]="$x2 $y";
   $hQrycontig{$tmpIdQrycontig}->[2]=$a[$ind_e1]; 
 
    #the matches line is red
   print O "2 1 0 2 4 0 40 0 -1 0.000 0 0 -1 0 0 2\n";
   print O "\t$x1 $y $x2 $y\n";
   print O "2 1 0 5 20 0 50 0 -1 0.000 0 0 -1 0 0 2\n";		
   printf O ("\t $x1 %.0f $x2 %.0f\n",$Yorig+150 , $Yorig+150);
}
#*********************************
sub print_connections{
 my ($setc1,$setx2,$sety2)=@_;                    # print "\nparam connect  @_\n";
 my ($setx1,$sety1) =split(/ /,$setc1);
 
   if ($Mgaps){      
     if ($setx1>$setx2){
        $tmpsetx1=$setx1;
        $setx1=$setx2;
        $setx2=$tmpsetx1;
     }
       $distx1x2=int(($setx2-$setx1)/2);
       $xcenter= $setx1+$distx1x2;
      
       if ($setx2-$setx1>4000)  { #if the distance is to big then heigh of the arc is set to 20 
         $heightArcUp = 20*$yscale;  
         $yoffcenter=int((($distx1x2**2)+$heightArcUp**2)*(1/(2*$heightArcUp)))-$heightArcUp ; 
       }
       else{
         $heightArcUp = int (0.447 * $distx1x2);#sectorul de cerc la 1/3 din raza.
         $yoffcenter=2*$heightArcUp;
       }
         print O "5 1 0 2 $linkcolor 0 50 0 -1 0.000 0 0 0 0 ";
         printf O ("%.3f %.3f $setx1 $sety1 $xcenter %.0f $setx2 $sety1 \n",$xcenter,$sety1+$yoffcenter,$sety1-$heightArcUp); 
  }
  else{
    print O "2 1 0 1 $linkcolor 0 50 0 -1 0.000 0 0 -1 0 0 2\n";
    print O "\t".$setc1." $setx2 $sety2\n";
  
  }
}
#*********************************
sub print_genes_mr{ 
 %hLastOnLevel=(); 
 @g=split(/\s+/,$hContig_genes{$a[-2]});
 for ($i=1;$i<=$#g;$i++){
 $kname=$g[$i];
   $tmpx2=0;
   $y=$Yorig-100-200*$geneinfo{$kname}->[0];
   #print id gena     
    $xid =$startdrawX+ int($Utr5End{$kname}/$Xscale);
    print_Id_genes() if ($printIDgenes);
   #     
   for ($l=2;$l<5;$l++){
     @c=split(";",$geneinfo{$kname}->[$l] ) ; 
      if (@c){ #print "de unde?@c\n" if (($l==2) or ($l==4));
      $colorend=$color{$l};
      if ($geneinfo{$kname}->[1] eq "-"){
         if ($l==2) { $colorend=$color{"4"}; }
         elsif ($l==4) {$colorend=$color{"2"};}
      }       
      for ($k=0;$k<=$#c;$k++){         
         @e=split (" ",$c[$k]); 
         $x1=$startdrawX+int($e[0]/$Xscale);
         $x2=$startdrawX+int($e[1]/$Xscale);
         if($x1==$x2) { $x2++;}
         if ( ($tmpx2) and ($x1-$tmpx2>1)){ #print the intron 
            print O "2 1 0 1 0 0 50 0 -1 0.000 0 0 -1 0 0 2\n";
            print O "\t $tmpx2 $y $x1 $y\n";
         }
         $tmpx2=$x2;
         print O "2 1 0 5 $colorend 0 50 0 -1 0.000 0 0 -1 0 0 2\n";#		
         print O "\t $x1 $y $x2 $y\n";  
      }
     }       
   }          
   #delete ($geneinfo{$kname});
 }
}

#*********************************
sub print_Id_genes{
 if (exists $hLastOnLevel{$geneinfo{$kname}->[0]}){
    $lastOnlevel=$hLastOnLevel{$geneinfo{$kname}->[0]};
    $printidspace = int(($Utr5End{$kname}-$Utr5End{$lastOnlevel})/$Xscale);
 }
 else{$printidspace=601;}
 if ($printidspace > 600){                
    #print contig name#
    print O "4 0 0 5 0 0 6 0.0000 4 90 270 ";
    printf O ("\t%.0f %.0f ",$xid,$y-50); 
    print O  $kname, "\\001\n";
    $hLastOnLevel{$geneinfo{$kname}->[0]}=$kname;
 } 
}
#*********************************
sub print_genes{
 %hLastOnLevel=();                                 
 foreach $kname  (sort {$Utr5End{$a} <=> $Utr5End{$b}} keys %Utr5End){  
   $tmpx2=0;
   if ($Utr5End{$kname}>$startcoord && $Utr5End{$kname}<$endcoord){       
     $y=$Yorig-100-200*$geneinfo{$kname}->[0];
     #print id gena     
        $xid = int(($Utr5End{$kname}-$startcoord)/$Xscale);
        print_Id_genes() if ($printIDgenes);
     #     
      for ($l=2;$l<5;$l++){
         @c=split(";",$geneinfo{$kname}->[$l] ); 
         $colorend=$color{$l};
         if ($geneinfo{$kname}->[1] eq "-"){
            if ($l==2) { $colorend=$color{"4"}; }
            elsif ($l==4) {$colorend=$color{"2"};}
         }
         
         for ($k=0;$k<=$#c;$k++){         
            @e=split (" ",$c[$k]); 
            $x1=int(($e[0]-$startcoord)/$Xscale);
            $x2=int(($e[1]-$startcoord)/$Xscale);
            if($x1==$x2) { $x2++;}
            if ( ($tmpx2) and ($x1-$tmpx2>1)){ #print the intron 
               print O "2 1 0 1 0 0 50 0 -1 0.000 0 0 -1 0 0 2\n";
               print O "\t $tmpx2 $y $x1 $y\n";
            }
            $tmpx2=$x2;
            print O "2 1 0 5 $colorend 0 50 0 -1 0.000 0 0 -1 0 0 2\n";#		
            print O "\t $x1 $y $x2 $y\n";  
         }      
      }
    }#endif "is in interval"
  }
 # delete ($Utr5End{$kname});
 # delete ($geneinfo{$kname});
 
}
#*********************************
sub print_grid{  

my ($xs,$xe,$startcontg,$endcontg)=@_;

$XlastPrint=0;$YlastPrint=0;

   #print ref contig
    print O "2 1 0 10 11 0 50 0 -1 0.000 0 0 -1 0 0 2\n";
    printf O ("\t $xs %.0f $xe %.0f\n",$Yorig+50,$Yorig+50);
    #print orizontal axes for PId (100%,75%,50%)
    for ($percent_id = 50; $percent_id < 101; $percent_id += 25) {
        print O  "2 1 2 1 0 7 60 0 -1 4.000 0 0 -1 0 0 2\n";
        printf O  ("\t$xs %.0f $xe %.0f\n",$Yorig+250+$YdistPID-($percent_id - 50) * $yscale,$Yorig+250+$YdistPID-($percent_id - 50) * $yscale);
        #last if ($Mgaps);    
    }
    #print orizontal markers for bp.
    $increment=10000/$Xscale;
    $no_incr=0;
    $xmark = $xs ;$xmark_float= $xs;
     while ($xmark < $xe){
        print O  "2 1 0 1 0 7 60 0 -1 0.000 0 0 -1 0 0 2\n";
        printf O  ("\t$xmark %.0f $xmark %.0f\n",$Yorig+$YdistPID+250,$Yorig+$YdistPID+300);
         #bp scale
        print O "4 0 0 100 0 0 8 0.0000 4 135 405 ";
        printf O ("\t %.0f %.0f",$xmark,$Yorig+$YdistPID+400); 
        print O  " $no_incr"."k", "\\001\n";
        $no_incr += 10;
        $xmark_float += $increment;
        $xmark=int($xmark_float);
    }
    
    #coord for chr ends
    print O "4 0 0 50 0 0 14 0.0000 4 135 450 $xs $Yorig $startcontg\\001\n";
    printf O ("4 0 0 50 0 0 14 0.0000 4 135 810 %.0f $Yorig $endcontg\\001\n",$xe-length($xe)*125);
     
    #print contig name#
    if ($mref){
      print O "4 0 0 5 0 0 8 0.0000 4 135 405 ";
      printf O ("\t%.0f %.0f ",$xs,$Yorig+70); 
      print O  $a[-2], "\\001\n";                     
    }
    #print vertical markers for PId scale 
    if (!$Mgaps){
      for ($percent_id = 50; $percent_id < 101; $percent_id += 25) {
        #left
        print O "4 0 0 100 0 0 8 0.0000 4 135 405 ";
        printf O ("\t%.0f %.0f", $xs-200,$Yorig+$YdistPID+250-($percent_id - 50) * $yscale + 20); 
        print O  " $percent_id%", "\\001\n";
        #right
        print O "4 0 0 100 0 0 8 0.0000 4 135 405 ";
        printf O ("\t%.0f %.0f",$xe+20, $Yorig+$YdistPID+250-($percent_id - 50) * $yscale+20 ); 
        print O  " $percent_id%", "\\001\n";
      
        # print the tick mark
        #left
       # print O  "2 1 0 1 0 7 60 0 -1 0.000 0 0 -1 0 0 2\n"; 
      #  printf O ("\t%.0f %.0f $xs %.0f\n",$xs-50,
	#   $Yorig+$YdistPID+250-($percent_id - 50) * $yscale, $Yorig+$YdistPID+250-($percent_id - 50) * $yscale);
      #right
       # print O  "2 1 0 1 0 7 60 0 -1 0.000 0 0 -1 0 0 2\n"; 
      #  printf O ("\t$xe %.0f %.0f %.0f\n",
	#   $Yorig+$YdistPID+250-($percent_id - 50) * $yscale,$xe+50, $Yorig+$YdistPID+250-($percent_id - 50) * $yscale);
      }
    }
    else{ # for Mgaps
      print O "4 0 0 100 0 0 7 1.5710 4 135 405 ";
      printf O ("\t%.0f %.0f", $xs-50,$Yorig+$YdistPID+250 - 5 * $yscale + 10); 
      print O  " + qry strand", "\\001\n";

      print O "4 0 0 100 0 0 7 1.5710 4 135 405 ";
      printf O ("\t%.0f %.0f", $xs-50,$Yorig+$YdistPID+250 - 30 * $yscale + 10); 
      print O  " - qry strand", "\\001\n";      
    }
    
}
#*********************************
sub print_legend{ 

     print O "4 0 0 100 0 0 8 0.0000 4 135 405 ";
     printf O ("\t%.0f %.0f ",100,$Yorig+$YdistPID+1100); 
     print O  " Legend ", "\\001\n";
     $y= $Yorig+$YdistPID+1300; #utr
     print O "2 1 0 1 0 0 50 0 -1 0.000 0 0 -1 0 0 2\n";#intron
     print O "\t 70 $y 99 $y\n";
     print O "2 1 0 5 27 0 50 0 -1 0.000 0 0 -1 0 0 2\n";
     print O "\t 100 $y 200 $y\n";
     print O "2 1 0 1 0 0 50 0 -1 0.000 0 0 -1 0 0 2\n";#intron
     print O "\t 200 $y 230 $y\n";        
     print O "4 0 0 100 0 0 8 0.0000 4 135 405 ";
     printf O ("\t%.0f %.0f ",300,$y+30); 
     print O  " 5' utr ", "\\001\n";
     $y += 150 ;#cds
     print O "2 1 0 1 0 0 50 0 -1 0.000 0 0 -1 0 0 2\n";#intron
     print O "\t 70 $y 99 $y\n";     
     print O "2 1 0 5 2 0 50 0 -1 0.000 0 0 -1 0 0 2\n";
     print O "\t 100 $y 200 $y\n";
     print O "2 1 0 1 0 0 50 0 -1 0.000 0 0 -1 0 0 2\n";#intron
     print O "\t 200 $y 230 $y\n";                 
     print O "4 0 0 100 0 0 8 0.0000 4 135 405 ";
     printf O ("\t%.0f %.0f ",300,$y+30); 
     print O  " cds ", "\\001\n";
     $y += 150;        #3' utr
     print O "2 1 0 1 0 0 50 0 -1 0.000 0 0 -1 0 0 2\n";#intron
     print O "\t 70 $y 99 $y\n";     
     print O "2 1 0 5 1 0 50 0 -1 0.000 0 0 -1 0 0 2\n";
     print O "\t 100 $y 200 $y\n";         
     print O "2 1 0 1 0 0 50 0 -1 0.000 0 0 -1 0 0 2\n";#intron
     print O "\t 200 $y 230 $y\n";             
     print O "4 0 0 100 0 0 8 0.0000 4 135 405 ";     
     printf O ("\t%.0f %.0f ",300,$y+30); 
     print O  " 3' utr ", "\\001\n";
     $y += 150; # match
     print O "2 1 0 2 4 0 50 0 -1 0.000 0 0 -1 0 0 2\n";
     print O "\t100 $y 200 $y\n";
     print O "4 0 0 100 0 0 8 0.0000 4 135 405 ";    
     printf O ("\t%.0f %.0f ",300,$y+30); 
     print O  " match found by $prog ", "\\001\n";
}
#*********************************
sub change_file_format{   

  $procfile =~ /^tmp(.+)/;
  $outfile = $1.$nrf.".".$format; 
  $comand = "fig2dev -L $format -x 100";
  $comand .= " -m $magn" if ($magn);
  $comand .= " -M ".$procfile.$nrf.".fig ".$outfile;
   
  $status =system($comand);
   print E "ERROR 1: fig2dev !\n" unless $status == 0;
   
  $status =system("rm $procfile".$nrf.".fig");
  
  if ($verb){
    print "$outfile";
    if ($mref){
      print "\n" ;
    }
    else {
      print "\t range : $startcoord\t$endcoord \n" ;
    }
  }
}
#*********************************
sub format_mgaps{
$tmpfile="tmpmgaps";
$tmpfile2=$alignm."coords" ;
get_ref_len();    #print $maxlenref."\n";
 open(M,">".$tmpfile2) || die "can't open \" $tmpfile2  \" file !";
  print M "$alignm\n";
  print M "Mgaps\n\n";
  print M "    [S1]     [E1]  |     [S2]     [E2]  |  [LEN 1]  [LEN 2]  |  [% IDY]  |  [LEN R]  [LEN Q]  |  [COV R]  [COV Q]  | [TAGS]\n";
  print M "===============================================================================================================================\n";

 open(T,">".$tmpfile) || die "can't open \" $tmpfile  \" file !";
 
 open(A,"<".$alignm) || die "can't open \" $alignm  \" file !";
 
 
 while(<A>) {
  chomp;
  @a=split;
  if ($a[0] =~ /^>/){    
    $nr_cluster=1;    
    $idquery=$a[1];
    if ($a[2] eq "Reverse"){$idquery .= "_rev";}
  }
  #elsif ($a[0] eq "#") {$nr_cluster++;}
  elsif($a[0] ne "#"){  
    $e1=$a[0]+$a[2];
    print T $a[0]."\t".$e1."\t"."|";
    print T "\t-\t-\t|"; 
    print T "\t-\t-\t|";
    print T "\t-\t|"; #pid
    print T "\t$maxlenref\t-\t|";#len seqs
   # print "\t-\t-\t|";
   # print "\t-\t-\t|";
   # print T " $Mgaps\t$idquery.$nr_cluster\n"; 
   print T " $Mgaps\t$idquery\n";   
  }
  
 }
 close(A);
 close(T);
 $command="sort -n -k 1 $tmpfile >> $tmpfile2";
 $status =system($command);  
 system("rm $tmpfile");
 close(M);
 $alignm=$tmpfile2;
 print STDERR "ERROR 1:  can't sort $tmpfile \n" unless $status == 0;
 print STDERR "\n**************************************** \n";
 print STDERR "New input file created : $alignm\n";
 print STDERR "**************************************** \n\n";
 
}
#*****************************
sub get_ref_len{ 
 $firstrow=1;
 open(A,"<".$alignm) || die "can't open \" $alignm  \" file !";
 while(<A>) {
   chomp;
   @a=split;
   if ($firstrow) {
     $firstrow=0;
     if ($a[0] !~ /^>/){
        print "\nWrong file format for MGAPS file : $alignm ! \n";
        exit;
     }
  }

  if ($a[0] =~ /^>/){ next; }
  elsif($a[0] ne "#"){  
    $e1=$a[0]+$a[2];
    $maxlenref=($maxlenref < $e1 ? $e1 : $maxlenref);
  }
 }
 close(A);
}
