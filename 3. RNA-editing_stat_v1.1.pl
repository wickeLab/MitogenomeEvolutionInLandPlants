use strict;
use warnings;
my $version="1.1 version";
use Getopt::Long;
use Text::CSV_PP;
my $parser = Text::CSV_PP->new();

###Written by: Yanlei Feng & Susann Wicke
###CITE: Feng & Wicke (2022) First mitochondrial genomes of true ferns allow modeling the mitogenomic inflation syndrome across all land plant lineages. bioarxiv.

my $usage = '
	Usage: perl this.pl -h -t table_file -f int -c int -o 
	   
	-t|--table: 	RNA editing table from GENEIOUS
		Column: all. Order doesn\'t matter
	-f|--frequency	frequency cutoff. 0-100, default: 0
	-c|--coverage	coverage cutoff. Default: 0
	-o|--other		also count u-to-c type
		if wants more type please change the code
';

my ( $help, $table, $other );
my $frequency = 0;
my $coverage = 0;
GetOptions( "help|h!" => \$help, 
			"table|t=s" => \$table,
			"frequency|f=s" => \$frequency,
			"coverage|c=s" => \$coverage,
			"other|o!" => \$other,
);
print "\n*************\n*$version*\n*************\n";

if( $help || ! -e $table ) {
	die $usage;
}

open IN, $table or die $!;
my $head = <IN>;
my %column;
$column{ $_ } = scalar keys %column for split /,/, $head;
my ( $gene_id, $species_id, $aa_change_id, $nt_change_id, $length_id, $substitution_id, $codon_position_id, $location_id, $coverage_id, $frequency_id );
if( exists $column{ 'CDS' } ){
	$gene_id = $column{ 'CDS' };	# atp1 CDS
}else{
	die "There's no CDS in $table!\n";
}
if( exists $column{ 'Length' } ){
	$length_id  = $column{ 'Length' }; # 1, 2, 3
}else{
	die "There's no Length in $table!\n";
}
if( exists $column{ 'Sequence Name' } ){
	$species_id = $column{ 'Sequence Name' }; # Orobanche_crenata-MT-1 - nad3 CDS - nad3 CDS
}else{
	die "There's no Sequence Name in $table!\n";
}
if( exists $column{ 'Amino Acid Change' } ){
	$aa_change_id = $column{ 'Amino Acid Change' }; # only when substitution, P -> S
}else{
	die "There's no Amino Acid Change in $table!\n";
}
if( exists $column{ 'Change' } ){
	$nt_change_id = $column{ 'Change' }; # C -> T, TAGG -> GTGC
}else{
	die "There's no Change in $table!\n";
}
if( exists $column{ 'Protein Effect' } ){
	$substitution_id = $column{ 'Protein Effect' }; # Substitution, Truncation, None and others (Extension, Start Codon Loss, Frame Shift)
}else{
	die "There's no Protein Effect in $table!\n";
}
if( exists $column{ 'CDS Position Within Codon' } ){
	$codon_position_id = $column{ 'CDS Position Within Codon' };
}else{
	die "There's no CDS Position Within Codon in $table!\n";
}
if( exists $column{ 'Coverage' } ){
	$coverage_id = $column{ 'Coverage' };
}else{
	die "There's no Coverage in $table!\n";
}
if( exists $column{ 'Minimum' } ){
	$location_id = $column{ 'Minimum' };
}else{
	die "There's no Minimum in $table!\n";
}
if( exists $column{ 'Variant Frequency' } ){
	$frequency_id = $column{ 'Variant Frequency' };
}else{
	die "There's no Variant Frequency in $table!\n";
}

foreach my $id ( $gene_id, $species_id, $aa_change_id, $nt_change_id, $length_id, $substitution_id, $codon_position_id, $location_id, $coverage_id, $frequency_id ){
	print "$id\n";
	die "Some parameter are problematic!\n" unless $id ne "";
}

print "Table looks fine. Processing...\n";

my %hash;
my %genes;
while( <IN> ){
	s/\r?\n//;
	$parser->parse( $_ );
	my @a = $parser->fields();
	
	$a[ $coverage_id ] =~ s/^(\d+).*/$1/;
	next if $a[ $coverage_id ] < $coverage;
	
	$a[ $frequency_id] =~ s/\%//g;
	$a[ $frequency_id] =~ s/^(\d+).*/$1/;
	next if $a[ $frequency_id] < $frequency;
	
	next unless $a[ $substitution_id ] =~ /Substitution|Truncation|None|Extension/;
	
	next unless $a[ $nt_change_id ] =~ /^[CT]+ -> [TC]+$/i;	# both U-to-C and C-to-U, if wants more type expression needs to be changed
	
	my $species = $a[ $species_id ] =~ s/([a-z]+[ _][a-z]+).*/$1/ri;
	my $gene = $a[ $gene_id ] =~ s/ CDS//ri;
	
	if( $a[ $length_id ] == 1 ){
		next if exists $hash{ $species }{ $gene }{ $a[ $location_id ] };
		my( $old, $new ) = $a[ $nt_change_id ] =~ /(\w) -> (\w)/;
		if( !$other ){
			next if $old eq 'T';
		}
		my $type = $old . "-" . $new;
		$hash{ $species }{ $gene }{ $type }{ 'total' }++;
		$hash{ $species }{ $gene }{ $type }{ 'pos1' }++ if $a[ $codon_position_id ] == 1;
		$hash{ $species }{ $gene }{ $type }{ 'pos2' }++ if $a[ $codon_position_id ] == 2;
		if( $a[ $substitution_id ] eq 'None' ){
			$hash{ $species }{ $gene }{ $type }{ 'silent' }++;
		}else{
			$hash{ $species }{ $gene }{ $type }{ 'nonsilent' }++;
		}
		$hash{ $species }{ $gene }{ $type }{ $a[ $location_id ] }++;
	}elsif( $a[ $length_id ] > 1 ){
		my( $old, $new ) = $a[ $nt_change_id ] =~ /^(\w+) -> (\w+)$/;
		next if length $old != length $new;
		my @old = split //, $old;
		my @new = split //, $new;
		my $maker = 1;
		for( my $i = 0; $i < scalar @old; $i++ ){
			last unless $old[ $i ];
			next if $old[ $i ] eq $new[ $i ];
			if( !$other ){
				next unless $old[ $i ] eq 'C';
				next unless $new[ $i ] eq 'T';
			}
			my $type = $old[ $i ] . "-" . $new[ $i ];
			next if exists $hash{ $species }{ $gene }{$type}{ $a[ $location_id ] + $i };
			$hash{ $species }{ $gene }{ $type }{ 'total' }++;
			$hash{ $species }{ $gene }{ $type }{ 'pos1' }++ if ( $a[ $location_id ] + $i ) % 3 == 1;
			$hash{ $species }{ $gene }{ $type }{ 'pos2' }++ if ( $a[ $location_id ] + $i ) % 3 == 2;
			
			if( $maker == 1 ){
				if( $a[ $substitution_id ] eq 'None' ){
					$hash{ $species }{ $gene }{ $type }{ 'silent' }++;
				}elsif( $a[ $substitution_id ] eq 'Substitution' ){
					$a[ $aa_change_id ] =~ s/ -> .*//;
					$hash{ $species }{ $gene }{ $type }{ 'nonsilent' } += length $a[ $aa_change_id ];
				}else{
					$hash{ $species }{ $gene }{ $type }{ 'nonsilent' }++;
				}
			}
			$hash{ $species }{ $gene }{ $type }{ $a[ $location_id ] + $i }++;
			$maker++;
		}
	}
	$hash{ $species }{ $gene }{ 'coverage_num' }++;
	$hash{ $species }{ $gene }{ 'coverage' } += $a[ $coverage_id];
	$genes{ $gene }++;
	
	# print "$species, $gene, $a[ $nt_change_id ], $a[ $substitution_id ], $a[ $codon_position_id ]\n";
}

my $output = $table =~ s/\.csv/\_output.csv/ri;
open OUT, ">$output";
foreach my $gene ( sort keys %genes ){
	print OUT ",$gene";
}
print OUT ",Sites,POS_1,POS_2,Synonymous,Non-Synonymous\n";

foreach my $species ( sort keys %hash ){
	print OUT "$species";
	my( $total_sites, $total_pos1, $total_pos2, $total_nonsilent, $total_silent );
	foreach my $gene ( sort keys %genes ){
		print OUT ",";
		if( exists $hash{ $species }{ $gene } ){
			my $cov = int( $hash{ $species }{ $gene }{ 'coverage' } / $hash{ $species }{ $gene }{ 'coverage_num' } );
			delete $hash{ $species }{ $gene }{ 'coverage' };
			delete $hash{ $species }{ $gene }{ 'coverage_num' };
			if( scalar keys %{ $hash{ $species }{ $gene } } > 0 ){
				foreach my $type ( sort keys %{ $hash{ $species }{ $gene } } ){
					$total_sites += $hash{ $species }{ $gene }{ $type }{ 'total' };
					$total_pos1 += $hash{ $species }{ $gene }{ $type }{ 'pos1' };
					$total_pos2 += $hash{ $species }{ $gene }{ $type }{ 'pos2' };
					$total_nonsilent += $hash{ $species }{ $gene }{ $type }{ 'silent' };
					$total_silent += $hash{ $species }{ $gene }{ $type }{ 'nonsilent' };
					print OUT "$type\{$hash{ $species }{ $gene }{ $type }{ 'total' }\($hash{ $species }{ $gene }{ $type }{ 'pos1' }\/$hash{ $species }{ $gene }{ $type }{ 'pos2' }\)\($hash{ $species }{ $gene }{ $type }{ 'nonsilent' }\/$hash{ $species }{ $gene }{ $type }{ 'silent' }\)\}";
				}
				print OUT "|$cov";
			}else{
				print OUT "0";
			}
		}else{
			print OUT "0";
		}
	}
	print OUT ",$total_sites,$total_pos1,$total_pos2,$total_nonsilent,$total_silent\n";
}

