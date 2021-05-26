use warnings;
use strict;

my ($ref, $isnvFile, $snpFile, $thres, $rs, $re, $idFile, $refname) = load_arg(\@ARGV);




if ($ref eq 'NULL') {help_info(); die};
if ($isnvFile eq 'NULL' && $snpFile eq 'NULL') {help_info(); die};
#if ($isnvFile eq 'NULL'){print "Warning: no isnv info file input, use only the SNP info.\n";} 
#if ($snpFile  eq 'NULL'){print "Warning: no snp info file input, use only the SNP info.\n";} 

my ($refHeader, $refSeq) = load_ref_fa($ref);

my $rh_mut = load_info($isnvFile, $snpFile, $thres);


# ------ print reference genome or genome region ----------
if ($refname eq 'NULL'){
    print ">$refHeader\n";
} else {
    print ">$refname\n";
}
my $refRegion = ref_region($refSeq, $rs, $re);
print "$refRegion\n";

# ------------------------------------
my $rh_id2newid;
$rh_id2newid = id2newid($idFile) if ($idFile ne 'NULL');



while (my ($samp, $rh_pos2mut) = each %$rh_mut){
    my $mutant = substitute_ref($refSeq, $rh_pos2mut, $rs, $re);    
    if ($idFile eq 'NULL'){
        print ">$samp\n";
    } else {
        if (defined $rh_id2newid->{$samp}){
            print ">$rh_id2newid->{$samp}\n";
        } else {
            die "[error] Cannot find the new id for $samp in $idFile\n";
        }
    }
    print "$mutant\n";
}


# ---------------------------------------------------------

sub help_info{
    print "\n\t**************************************************************************************\n";
    print "\t* Generate the consenus sequences by integrating reference geonome and iSNV&SNP info *\n";
    print "\t* from step5 of iSNV_calling pipeline.                                               *\n";
    print "\t*                                                                       nim, 2017-10 *\n";
    print "\t**************************************************************************************\n";
    
    print "\nUsage:\n  perl genConsensusFa.pl -r <ref.fa> -isnv <iSNV_info.txt> -snp <SNP_info.txt> [-thres freq_threshold -rs <start> -re <end> -id <id2newid file> -refname <name>] \n";
     print "\n  Options:\n";
     print "  \t-r           fasta file of the reference genome.\n";
     print "  \t-isnv        isnv_info.txt output by the step5 script of iSNV_calling pipeline.\n";
     print "  \t-isnp        isnp_info.txt output by the step5 script of iSNV_calling pipeline.\n";
     print "  \t-thres       The threshold of the isnv freq for substitution. Default = 0.5\n";
     print "  \t-rs          The start location of reference to generate consensus sequences.\n";
     print "  \t-re          The end location of reference to generate consensus sequences.\n";
     print "  \t-id          The file contains the id to new id, format: ID\tNewID.\n";
     print "  \t-refname     User customed reference name. Default = the header of reference fasta.\n";
}

sub load_arg{
    my $ra      = shift;
    my $ref     = 'NULL';
    my $isnv    = 'NULL';
    my $snp     = 'NULL';
    my $rs      = 'NULL';
    my $re      = 'NULL';
    my $thres   = 0.5;
    my $idfile  = 'NULL';
    my $refname = 'NULL';
    for(my $i = 0; $i < @{$ra} - 1; $i+=2){
        my $ar = $$ra[$i];
        if ($ar eq '-r'){
            $ref = $$ra[$i+1];
        } elsif ($ar eq '-isnv'){
            $isnv = $$ra[$i+1];
        } elsif ($ar eq '-snp'){
            $snp = $$ra[$i+1];
        } elsif ($ar eq '-rs'){
            $rs = $$ra[$i+1];
        } elsif ($ar eq '-re'){
            $re = $$ra[$i+1];
        } elsif ($ar eq '-thres'){
            $thres = $$ra[$i+1];
        } elsif ($ar eq '-id'){
            $idfile = $$ra[$i+1];
        } elsif ($ar eq '-thres'){
            $thres = $$ra[$i+1]; 
        } elsif ($ar eq '-refname'){
            $refname = $$ra[$i+1];
        } else {
            help_info();
            die "\n[error] no option for \"$ar\".\n";
        }
    }
    return ($ref, $isnv, $snp, $thres, $rs, $re, $idfile, $refname);
}

sub load_info{
    my $file1 = shift;
    my $file2 = shift;
    my $t = shift;
    my $rh;
    if ($file1 ne 'NULL'){
        open FP, "$file1" or die "$!: $file1\n   ";
        while (<FP>){
            chomp;
            next if (/^#/ || length($_) < 50);
            my @row = split /\s+/, $_;
            my $samp= $row[0];
            my $pos = $row[1];
            my $MuAF= $row[3];
            $row[4] =~ /(\w):(\w)-(\w)/;
            my $ref = $1; 
            my $max = $2;
            my $sec = $3;
            my $mut = $max eq $ref? $sec : $max;
            if ($MuAF >= $t){
                $rh->{$samp}->{$pos} = $mut;
            }     
        }
        close FP;
    }
    if ($file2 ne 'NULL'){
        open FP, "$file2" or die "$!: $file2\n   ";
        while (<FP>){
            chomp;
            next if (/^#/ || length($_) < 50);
            my @row = split /\s+/, $_;
            my $samp= $row[0];
            my $pos = $row[1];
            my $MuAF= $row[3];
            $row[4] =~ /(\w):(\w)-(\w)/;
            my $ref = $1; 
            my $max = $2;
            my $sec = $3;
            my $mut = $max eq $ref? $sec : $max;
            if ($MuAF >= $t){
                $rh->{$samp}->{$pos} = $mut;
            }     
        }
        close FP;
    } 
    return $rh;
}

sub load_ref_fa{
    open FP, "$_[0]" or die "$!: $_[0]\n   ";
    my $header;
    my $seq = '';
    while (<FP>){
        chomp;
        if(/^>(.*)/){
            $header = $1;
        } elsif (/^\w+$/){
            $seq .= $_;
        }
    }
    close FP;
    return ($header, $seq);
}

sub substitute_ref{
    my $ref = shift;
    my $rh = shift;
    my $s = shift;
    my $e = shift;
    my $mutant = '';
    my @mutant = split //, $ref;
    my @pos = sort {$a<=>$b} keys %$rh;
    for(my $i = 0; $i < @pos; $i++){
        $mutant[$pos[$i]-1] = $rh->{$pos[$i]};
    }
    my $totbp = @mutant;
    if ($s ne 'NULL' && $e ne 'NULL'){
        @mutant = @mutant[$s-1..$e-1];
    } elsif ($s ne 'NULL') {
        @mutant = @mutant[$s..$totbp-1];
    } elsif ($e ne 'NULL') {
        @mutant = @mutant[0..$e];
    }
    return join('',@mutant);
}


sub ref_region{
    my $seq = shift;
    my $s = shift;
    my $e = shift;
    my $reg;
    my $totbp = length($seq);
    if ($s ne 'NULL' && $e ne 'NULL'){
        $reg  = substr($seq, $s-1, $e-$s+1);
    } elsif ($s ne 'NULL') {
        $reg  = substr($seq, $s-1);
    } elsif ($e ne 'NULL') {
        $reg  = substr($seq, 0, $e-1);
    } else {
        return $seq;
    }
    return $reg;
}

sub id2newid{
    open FP, "$_[0]" or die "$!: $_[0]\n";
    my $rh;
    while (<FP>){
        chomp;
        my @row = split /\s+/, $_;
        next if (@row != 2);
        $rh->{$row[0]} = $row[1];
    }
    close FP;
    return $rh;
}

