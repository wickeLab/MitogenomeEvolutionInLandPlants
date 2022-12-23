use strict;
use warnings;
use File::Basename;

#### this script reads blastn table (mitogenome against pt database) and counts pt-derived length and number, also gives the annotation table of pt-derived.

###Written by: Yanlei Feng & Susann Wicke
###CITE: Feng & Wicke (2022) First mitochondrial genomes of true ferns allow modeling the mitogenomic inflation syndrome across all land plant lineages. bioarxiv.

my $length_cutoff = 100; # base pair
my $identity_cutoff = 0; # %

if( $ARGV[0] && $ARGV[0] =~ /^\d+$/ ){
	$length_cutoff = $ARGV[0];
}

print "repeat counting: use pt-derived length threshold: $length_cutoff\n";

my $out_table = "plastid_derived_stat.csv";
open STAT, ">$out_table" or die $!;
print STAT "Species,Total pt length,pt number\n";

my @xls=<xls/*.xls>;			# find the BLASTN tables under the folder
foreach my $file ( @xls ){
	my $base = basename $file;
	my $species = $base =~ s/\.fasta//r;
	print $species, "\n";
	
	my %NumPair;
	open TA, $file or die $!;
	while( <TA> ){
		s/\r?\n//;
		my @a = split /\t+/;
		
		next if $a[3] < $length_cutoff;			# skip hits shorter than the length cutoff
		next if $a[2] < $identity_cutoff;		# skip identity less than cutoff
				
		for( my $j = $a[6]; $j <= $a[7]; $j++ ){
			$NumPair{ $a[0] }{ $j }++;
		}
	}
	close TA;
	
	# my $sequin = $base =~ s/\.xls/\.plastid/r;
	# open OUT, ">$sequin" or die $!;
	my $pt_num = 0; my $total_length;
	my $i = 1;
	foreach my $key ( sort keys %NumPair ){
		
		$total_length += scalar keys %{ $NumPair{ $key } };
		
		my %num = combin_num( keys %{ $NumPair{ $key } } );
		$pt_num += scalar keys %num;
		# my $j = 1;
		# foreach my $begin ( sort keys %num ){
		# 	if( $j == 1 ){
		# 		print OUT ">Features $key\n";
		# 		$j++;
		# 	}
		# 	print OUT "$begin\t$num{ $begin }\tmisc_feature\n";
		# 	print OUT "\t\t\tnote\tplastid-derived $i\n";
		# 	$i++;
		# }
	}
	# close OUT;
	print STAT "$species,$total_length,$pt_num\n";
}
close STAT;

sub combin_num{
	my @list=@_;
	@list=sort {$a<=>$b} @list;
	my $list=join(",",@list);
	1 while $list =~ s/(-)?(\d+),(\d++)(?(?{$3 != $2+1})(*F))/ defined $1 ? "-$3" : "$2-$3" /e;
	# return $list;
	my @get=split /,/,$list;
	@get=grep{/\d+\-\d+/}@get;
	my %out;
	foreach my $sp(@get){
		my($n1,$n2)=split/-/,$sp;
		$out{$n1}=$n2;
	}
	return %out;
}

