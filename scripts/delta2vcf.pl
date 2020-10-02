#!@PERL@
#this code converts delta file (stdin) to vcf file (stdout)
if($#ARGV>=0){
  die("delta to vcf converter. vcf file contains variants in reference detected by query alignments. The best practice is to ude delta-filter -1 before running delta2vcf, i.e. delta-filter -1 file.delta | delta2vcf > file.vcf.\nUsage: delta2vcf < file.delta > file.vcf\n");
}
$line=<STDIN>;
chomp($line);
my ($ref,$qry)=split(/\s+/,$line);
$line=<STDIN>;
chomp($line);
die ("not a valid delta file") unless($line eq "NUCMER");

#loading reference and querysequences
open(FILE,$ref);
while($line=<FILE>){
  chomp($line);
  if($line=~/^\>/){
    my @f=split(/\s+/,$line);
    if(not($seq eq "")){
      $rseq{$ctg}=$seq;
      $seq=reverse($seq);
      $seq=~tr/ACGTNacgtn/TGCANtgcan/;
      $rseq_rc{$ctg}=$seq;
    }
    $ctg=substr($f[0],1);
    $seq="";
  }else{
    $seq.=$line;
  }
}
if(not($seq eq "")){
  $rseq{$ctg}=$seq;
  $seq=reverse($seq);
  $seq=~tr/ACGTNacgtn/TGCANtgcan/;
  $rseq_rc{$ctg}=$seq;
}
open(FILE,$qry);
while($line=<FILE>){
  chomp($line);
  if($line=~/^\>/){
    my @f=split(/\s+/,$line);
    if(not($seq eq "")){
      $qseq{$ctg}=$seq;
      $seq=reverse($seq);
      $seq=~tr/ACGTNacgtn/TGCANtgcan/;
      $qseq_rc{$ctg}=$seq;
    }
    $ctg=substr($f[0],1);
    $seq="";
  }else{
    $seq.=$line;
  }
}
if(not($seq eq "")){
  $qseq{$ctg}=$seq;
  $seq=reverse($seq);
  $seq=~tr/ACGTNacgtn/TGCANtgcan/;
  $qseq_rc{$ctg}=$seq;
}
#print header
($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime();
print "##fileformat=VCFv4.2\n";
printf("##fileDate=20%02d%02d%02d\n",$year-100,$mon+1,$mday);
print "##source=delta2vcf\n";
print "##reference=file:/$ref\n";
foreach $c(keys %rseq){
  print "##contig=<ID=$c,length=",length($rseq{$c}),">\n";
}
print "##FORMAT=<ID=GT,Number=1,Type=Integer,Description=\"Genotype\">\n";
print "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Coverage Depth, 1 by default\">\n";
print "##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Number of observations for each allele, 1,2 by default\">\n";
print "##FORMAT=<ID=RO,Number=1,Type=Integer,Description=\"Reference allele observation count, 1 by default\">\n";
print "##FORMAT=<ID=QR,Number=1,Type=Integer,Description=\"Sum of quality of the reference observations, 40 by default\">\n";
print "##FORMAT=<ID=AO,Number=A,Type=Integer,Description=\"Alternate allele observation count, 2 by default\">\n";
print "##FORMAT=<ID=QA,Number=A,Type=Integer,Description=\"Sum of quality of the alternate observations, 40 by default\">\n";
print "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tCOUNTS\n";
@outarray=();
#going through the delta file
while($line=<STDIN>){
  chomp($line);
  my @f=split(/\s+/,$line);
  print "$line $#f\n" if($DEBUG);
  if($#f==3){#we found delta header
    $refname=substr($f[0],1);
    $qryname=$f[1];
    $reflen=$f[2];
    $qrylen=$f[3];
  }elsif($#f==6){#we found new alignment header
    $ref_start=$f[0];
    $ref_end=$f[1];
    $qry_start=$f[2];
    $qry_end=$f[3];
    @delta_begins=();
    @delta_lengths=();
    $last_indel=0;
    while($line=<STDIN>){
      print $line if($DEBUG);
      chomp($line);
      last if($line eq "0");
      if($line>1||$line<-1){
        if($line>0){
          push(@delta_begins,$line+$last_indel);
          push(@delta_lengths,1);
          $last_indel+=$line;
        }else{
          push(@delta_begins,-(-$line+$last_indel));
          push(@delta_lengths,1);
          $last_indel+=-$line;
        }
      }else{
        $delta_lengths[-1]++;
        $last_indel++;
      }
    }
#we now read the deltas for this alignment into the arrays.  Next line is either another delta or a new alignment
#let's process the deltas
    print "Found an alignment $refname $qryname $ref_start $ref_end $qry_start $qry_end\n" if($DEBUG);
    print join(" ",@delta_begins),"\n" if($DEBUG);
    print join(" ",@delta_lengths),"\n" if($DEBUG);
    $delta_index=0;
    $rseq=$rseq{$refname};
    if($qry_start<$qry_end){
      $qseq=$qseq{$qryname};
    }else{
      $qseq=$qseq_rc{$qryname};
      $qry_start=length($qseq)-$qry_start+1;
      $qry_end=length($qseq)-$qry_end-1;
    }
    $ref_pos=$ref_start;
    $qry_pos=$qry_start;
    $alignment_coord=1;
    while($ref_pos<=$ref_end){
      if($alignment_coord==abs($delta_begins[$delta_index])){#insertion or deletion
        if($delta_begins[$delta_index]>0){
          for($j=0;$j<$delta_lengths[$delta_index];$j++){
            print substr($rseq,$ref_pos-1,1)," *"," $ref_pos $qry_pos $alignment_coord\n" if($DEBUG);
            $ref_pos++;
            $alignment_coord++;
          }
          push @outarray, {'refname'=>$refname, 'refpos'=>($ref_pos-$delta_lengths[$delta_index]-1),'ref'=>substr($rseq,$ref_pos-$delta_lengths[$delta_index]-2,$delta_lengths[$delta_index]+1),'qry'=>substr($qseq,$qry_pos-2,1) };
          $delta_index++;
        }else{#deletion
          for($j=0;$j<$delta_lengths[$delta_index];$j++){
            print "* ",substr($qseq,$qry_pos-1,1)," $ref_pos $qry_pos $alignment_coord\n" if($DEBUG);
            $qry_pos++;
            $alignment_coord++;
          }
          push @outarray,{'refname'=>$refname,'refpos'=>($ref_pos-1),'ref'=>substr($rseq,$ref_pos-2,1),'qry'=>substr($qseq,$qry_pos-$delta_lengths[$delta_index]-2,$delta_lengths[$delta_index]+1)};
          $delta_index++;
        }
      }else{ 
        print substr($rseq,$ref_pos-1,1)," ",substr($qseq,$qry_pos-1,1)," $ref_pos $qry_pos $alignment_coord\n" if($DEBUG);
        push @outarray, {'refname'=>$refname,'refpos'=>$ref_pos,'ref'=>substr($rseq,$ref_pos-1,1),'qry'=>substr($qseq,$qry_pos-1,1)} if(not(uc(substr($rseq,$ref_pos-1,1)) eq uc(substr($qseq,$qry_pos-1,1))));
        $ref_pos++;
        $qry_pos++;
        $alignment_coord++;
      }
    }#while
  }#alignment
}#top loop

#sort output
@outarray_sorted=sort { $a->{'refname'} cmp $b->{'refname'} || $a->{'refpos'} <=> $b->{'refpos'} } (@outarray);
#print
$last_refname="";
$last_pos="";
foreach $l(@outarray_sorted){
print $l->{'refname'},"\t",$l->{'refpos'},"\t\.\t",$l->{'ref'},"\t",$l->{'qry'},"\t40\tPASS\t.\tGT:DP:AD:RO:QR:AO:QA\t0:1:1,2:1:40:2:40\n" if(not($last_refname eq $l->{'refname'}) || not($last_pos==$l->{'refpos'}));
$last_refname=$l->{'refname'};
$last_pos=$l->{'refpos'};
}

