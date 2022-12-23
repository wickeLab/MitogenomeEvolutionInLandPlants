# use strict;
use warnings;

#### this script count the real and total also number of duplicated part of the genome.

###Written by: Yanlei Feng & Susann Wicke
###CITE: Feng & Wicke (2022) First mitochondrial genomes of true ferns allow modeling the mitogenomic inflation syndrome across all land plant lineages. bioarxiv.


my $length_cutoff = 100; # base pair
my $identity_cutoff = 95; # %

if( $ARGV[0] && $ARGV[0] =~ /^\d+$/ ){
	$length_cutoff = $ARGV[0];
}

print "repeat counting: use repeat length threshold: $length_cutoff\n";

my @xls=<fasta/*.xls>;			# find the BLASTN tables under the folder
open OUT, ">repeat_stat.csv" or die $!;
print OUT "Species,Length,Repeat_Num_", $length_cutoff, ",Long_repeat (> 500B),Long_repeat (> 1K),Long_repeat (> 4K),Real\n";
foreach $file ( @xls ){
	my %count100;
	my %repeat100;
	my $r100 = 0;
	my $total = 0;
	my $gc_num = 0;
	my $long_4000 = 0;
	my $long_1000 = 0;
	my $long_500 = 0;
	my $long_100 = 0;
	my $gc_ratio;
	open IN, $file or die $!;
	my $species = $file;
	$species =~ s/\.\w+$//;
	$species =~ s/[_\-]annotated//;
	$species =~ s/-MT//;
	print "$species\n";
	
	### read the fasta file to get the sequence length
	my $fasta = $file =~ s/\.xls/\.fasta/r;
	open FA, $fasta or die $!;
	my %fasta;
	while( <FA> ){
		chomp( my $name = $_ );
		chomp( my $seq = <FA> );
		$name =~ s/>//;
		$name =~ s/\r//g;
		$seq =~ s/\r//g;
		$name =~ s/\s+.*//;
		$fasta{ $name } = $seq;
		$total += length $seq;
	}
	
	my %real;
	my $real_length = 0;
	
	while(<IN>){
		my @a = split /\t/;
		my $con = $a[0];
		next if $a[3] < $length_cutoff;			# skip hits shorter than the length 
		next if $a[3] eq length $fasta{ $con };	# skip complete match (themselves)
		next if $a[2] < $identity_cutoff;		# skip hit identity less than cutoff
		
		# count total length
		$r100++;
		$count100{ $con }++;
		for( my $j = $a[6]; $j <= $a[7]; $j++ ){
			$repeat100{ $con }{ $j }++;
		}
		
		if($a[3] > 999){
			$long_1000++;
		}
		if($a[3] > 499){
			$long_500++;
		}
		
		if($a[3] > 3999){
			$long_4000++;
		}
		
		# count real duplicated
		my $step1 = 0;
		foreach my $chrm ( sort keys %real ){
			$step1 += scalar keys %{ $real{ $chrm } };
		}
		 
		for( my $m = $a[6]; $m <= $a[7]; $m++ ){
			$real{ $con }{ $m }++;
		}
		
		my $step2 = 0;
		foreach my $chrm ( sort keys %real ){
			$step2 += scalar keys %{ $real{ $chrm } };
		}
		my $change1 = $step2 - $step1;
		
		if( $a[8] > $a[9] ){
			( $a[8], $a[9] ) = ( $a[9], $a[8] );
		}
		for( my $m = $a[8]; $m <= $a[9]; $m++ ){
			$real{ $a[1] }{ $m }++;
		}
		
		my $step3;
		foreach my $chrm ( sort keys %real ){
			$step3 += scalar keys %{ $real{ $chrm } };
		}
		my $change2 = $step3 - $step2;
		
		if( $change1 >= $change2 ){
			$real_length += $change1;
		}else{
			$real_length += $change2;
		}
	}
	
	# my $sum100;
	# foreach my $key ( sort keys %repeat100 ){
	# 	$sum100 += scalar keys %{ $repeat100{ $key } };
	# 	#my %num = combin_num( keys %{ $repeat100{ $key } } );
	# 	#foreach my $begin ( sort keys %num ){
	# 	#	my $sub = substr( $fasta{ $key }, $begin-1, $num{ $begin } - $begin + 1 );
	# 	#	my $gc = $sub =~ tr/GC/GC/;
	# 	#	$gc_num += $gc;
	# 	#}
	# }
	
	$r100 = int( $r100 / 2 );
	$long_500 = int( $long_500 / 2 );
	$long_1000 = int( $long_1000 / 2 );
	$long_4000 = int( $long_4000 / 2 );
	
	print OUT "$species,$total,$r100,$long_500,$long_1000,$long_4000,$real_length\n";
}

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
		# if($n2-$n1>=$threshold){		# change threshold
			$out{$n1}=$n2;
		# }
	}
	return %out;
}

