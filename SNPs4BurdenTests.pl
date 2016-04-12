use strict;
use warnings;
use Data::Dumper; # used for diagnostic purposes.
use Pod::Usage qw(pod2usage); # Used for providing useful error/warning and user messages.

# GENCODE: gene, exon, transcript, CDS, UTR
# GTEx: promoter, CTCF, enhancer, promoterFlank, openChrom, TF_bind, allreg
# overlap: promoter, CTCF, enhancer, promoterFlank, openChrom, TF_bind, allreg

# As a first step the name of the gene will be parsed to using the Ensembl Perl API to get the genomic location and the stable ID.
my $modules = join(" ", @INC);
unless ($modules =~ // and $modules =~ //) {
    print STDERR "[Warning] Ensembl packages were not available.\n[Warning] Required packages were loaded from /nfs/team144/software/ensembl-releases/83\n";
    use lib "/nfs/team144/software/ensembl-releases/83/bioperl-1.6.1";
    use lib "/nfs/team144/software/ensembl-releases/83/ensembl/modules";
    use lib "/nfs/team144/software/ensembl-releases/83/ensembl-variation/modules";
}

# Loading required modules:
use Bio::EnsEMBL::Registry;

# Initializing registry:
my $registry = 'Bio::EnsEMBL::Registry';
$registry->load_registry_from_db(
  -host => 'ensembldb.ensembl.org',
  -user => 'anonymous',
);

# Establishing connection with the corresponding databases:
our $ga = $registry->get_adaptor( "human", "core", "gene" );

# Files for each required dataset are hardcoded:
our $geneBedFile      = '/lustre/scratch113/projects/helic/ds26/project_burden/Annotation_file/Annotations_GENCODE_GTEx_overlap.bed.gz';
our $vcfFile          = '/lustre/scratch115/projects/t144_helic_15x/manolis-initial-385/15x_helic_manolis_initial_385_batch.vcf.gz';
our $chainFile        = '/nfs/team144/ds26/projects/2016.03.21_Regulation_GTEx/hg38ToHg19.over.chain.gz';
our $temporaryBedFile = '/nfs/team144/ds26/projects/2016.03.21_Regulation_GTEx/temporary_hits_GRCh38_%s.bed';
our $liftoverPath     = '/nfs/users/nfs_d/ds26/bin/liftOver';
our $EigenPath        = '/lustre/scratch113/teams/zeggini/users/ds26/refseq/eigen_v1.0/Eigen_hg19_0916_%s.tab.bgz';
our $caddPath         = '/lustre/scratch114/resources/cadd_scores/20150729-v1.3/whole_genome_SNVs_inclAnno.tsv.gz';

# All files must exist otherwise the scrip quits:
foreach my $file ($geneBedFile, $vcfFile,  $chainFile, $liftoverPath, $EigenPath, $caddPath){
    $file =~ s/\%s/chr12/;
    pod2usage({-verbose => 99, -message => "[Error] One of the requested file does not exits ($file). Exiting.\n", -sections => "FILES|SYNOPSIS" })unless ( -e $file);
}

# Initializing command line parameters:
use Getopt::Long qw(GetOptions);
my ($inputFile, $outputFile, $GENCODE, $GTEx, $verbose, $overlap, $help,
    $extend, $minor, $tissue_list, $MAF, $MAC, $score, $MAF_weight, $k);

GetOptions(
    # Input/Output:
    'input|i=s' => \$inputFile,
    'output|o=s' => \$outputFile,

    # GENCODE features:
    "GENCODE|G=s" => \$GENCODE,

    # Extend regions with a defined length:
    "extend=s" => \$extend,

    # Appris:
    "SkipMinor" => \$minor,

    # GTEx:
    "GTEx|E=s" => \$GTEx,
    "overlap|L=s" => \$overlap,

    # Only considering regulatory elements associated in the
    # following tissues:
    "tissues=s" => \$tissue_list,

    # Variant features:
    'MAF=s' => \$MAF,
    'MAC=s' => \$MAC,

    # Verbose output:
    'verbose|v' => \$verbose,

    # Which score we need:
    'score|s=s' => \$score,

    # A flag to weight by MAF:
    'MAFweight|w' => \$MAF_weight,
    'ExpConst|k=s' => \$k,

    # Asking for help:
    'help|h' => \$help,
);

# Checking essential parameters:
pod2usage({-verbose => 2}) if $help;
pod2usage({-message => "[Error] Input file has to be specified with the -i switch!\n", -verbose => 0}) unless $inputFile;
pod2usage({-message => "[Error] Output file prefix has to be specified with the -o switch!\n", -verbose => 0}) unless $outputFile;
pod2usage({-message => "[Warning] Scoring method was not specified! Eigen scores are used.\n", -verbose => 0, -exitval => "NOEXIT"}) unless $score;
pod2usage({-message => '[Warning] Minor allele frequency has to be between 0 and 1.\n', -exitval => "NOEXIT"}) unless $MAF > 0 and $MAF <= 1;

# Parsing comman line arguments:
our %GENCODE = %{&parseGENCODE($GENCODE)};
our %GTEx = %{&parseRegulation($GTEx)};
our %overlap = %{&parseRegulation($overlap)};
our %Variant = ('MAF' => $MAF || 0.05,
                'MAC' => $MAC || 0);

# Optional parameters regarding the genes:
$GENCODE{'minor'}  = 1 if $minor;
$GENCODE{'extend'} = $extend || 0;

# Assigning default value to the weighting constant:
$k = 50 unless $k;
$score = "Eigen" unless $score;
pod2usage({-message => sprintf('[Warning] %s is not a supported scoring method. Either CADD, GERP or Eigen has to be specified.
[Warning] The default Eigen method is used.\n', $score), -verbose => 0, -exitval => "NOEXIT"}) unless ( $score eq "GERP" or
                    $score eq "CADD" or $score eq "Eigen"); # The submitted scoring method might not be supported! Check for it now!


###
### Reading file, processing the list of genes, line by line.
###
open(my $INPUT, "<", $inputFile) or die "[Error] Input file ($inputFile) could not be opened.\n";
open(my $variant_output, ">", $outputFile."_variants") or die "[Error] Output file ($outputFile\_variants) could not opened.\n";
my $genotypes= {};
while (my $ID = <$INPUT>) {
    chomp $ID;
    my ($chr, $start, $end, $stable_ID, $name) = &GetCoordinates($ID);

    print "\n\n[Info] Queried gene: $name (Ensembl ID: $stable_ID), Genomic location: chr$chr:$start-$end\n";
    my $bedlines = &BedToolsQuery($chr, $start, $end, $stable_ID);

    my @filteredLines = &FilterLines($bedlines);
    print "[Info] Number of selected genomic regions: ", scalar(@filteredLines),"\n" if $verbose;

    # Collapsing bed string:
    my $CollapsedBed = &CollapseRegions(@filteredLines);

    # Retrieve overlapping variations:
    my $variants = &GetVariants($CollapsedBed);

    # Process variant list:
    my $hash;
    ($hash, $genotypes) = &processVar($variants, $genotypes);

    # Saving temporary bedfile for liftover:
    $hash = &liftover($hash, $name);

    # Returning Eigen scoresL:
    $hash = &get_Eigen_Score($hash) if $score eq "Eigen";

    # Get CADD and GERP scores:
    $hash = &get_CADD_GERP($hash, $score) if $score eq "GERP" or $score eq "CADD";

    # Weighting scores if requested:
    $hash = &weight_score($hash, $k) if $MAF_weight;

    # Saving data into file:
    &print_SNPlist($hash, $variant_output, $name);

    # Diagnostic dumper of the hash:
    # print Dumper $hash;
}

# print "Genotypes:", Dumper ($genotypes);
open(my $genotype_output, ">", $outputFile."_genotype") or die "[Error] Output file ($outputFile\_genotype) could not opened.\n";
&print_genotypes($genotypes, $genotype_output);

##
## Function to parsing input parameters
##
sub parseGENCODE {
    # Accepted features: gene, exon, transcript, CDS, UTR
    my $gencodeString = $_[0];
    my %hash = ();
    foreach my $feature (split(",", $gencodeString)){
        $hash{$feature} = 1;
    }
    return \%hash
}
sub parseRegulation {
    # Accepted features: promoter, CTCF, enhancer, promoterFlank, openChrom, TF_bind, allreg
    my $regString = $_[0];
    my %hash = ();
    foreach my $feature (split(",", $regString)){
        $hash{'TF binding site'}   = 1 if "TF_bind" eq $feature;
        $hash{'CTCF Binding Site'} = 1 if "CTCF" eq $feature;
        $hash{'Enhancer'}          = 1 if "enhancer" eq $feature;
        $hash{'Open chromatin'}    = 1 if "openChrom" eq $feature;
        $hash{'Promoter'}          = 1 if "promoter" eq $feature;
        $hash{'allreg'}            = 1 if "allreg" eq $feature;
        $hash{'Promoter Flanking Region'} = 1 if "promoterFlank" eq $feature;
    }
    return \%hash
}

##
## Saving genotype information:
##
sub print_genotypes  {
    my %genotype = %{$_[0]};
    my $outputhandler = $_[1];

    # Get list of sample IDs:
    my $samples = `zgrep -m1 CHROM  $vcfFile | cut -f10-`;

    # Saving data:
    print $outputhandler "0\t$samples";
    for my $var (keys %genotype){
        print $outputhandler "$var\t", join("\t", @{$genotype{$var}}), "\n";
    }
}

##
## Weight scores with MAF. (lower MAF higher weight.)
##
sub weight_score {
    my %hash = %{$_[0]};
    my $k = $_[1];

    foreach my $snpID (keys %hash){
        my $AC = $hash{$snpID}{'frequencies'}[0];
        my $AN = $hash{$snpID}{'frequencies'}[1];
        my $MAF = $hash{$snpID}{'frequencies'}[2];
        $hash{$snpID}{score} = $hash{$snpID}{score} * exp( - ($AC - 1)/$AN * $k);
    }
    print "[Info]";
    return \%hash;
}

##
## Saving variants in SNP list:
##
sub print_SNPlist {
    my %hash = %{$_[0]};
    my $outputhandle = $_[1];
    my $gene_name = $_[2];

    #
    print $outputhandle "$gene_name\t1\t", join("\t", keys %hash),"\n$gene_name\t0\t";
    foreach my $snpid (keys %hash){
        print  $outputhandle "$hash{$snpid}{score}\t"
    }
    print $outputhandle "\n";

    return 1;
}


##
## Get Eigen scores
##
sub get_Eigen_Score {
    my %hash = %{$_[0]};

    # Looping throuh all variants and return the Eigen score for all:
    foreach my $var (keys %hash){
        my $EigenFile = sprintf($EigenPath, $hash{$var}{GRCh37}[0]);
        (my $chr = $hash{$var}{GRCh37}[0] ) =~ s/chr//i;
        my $tabix_query = sprintf("tabix %s %s:%s-%s | grep %s", $EigenFile, $chr, $hash{$var}{GRCh37}[2], $hash{$var}{GRCh37}[2], $hash{$var}{alleles}[1]);
        my $lines = `bash -O extglob -c \'$tabix_query\'`;
        $hash{$var}{score} = "NA"; # Initialize Eigen score.

        foreach my $line (split("\n", $lines)){
            # Line: 12	56482614	C	A	1.42750841783102
            chomp $line;
            my ($chr, $pos, $ref, $alt, $score) = split("\t", $line);

            # Testing alleles:
            if (($ref eq $hash{$var}{alleles}[0] and $alt eq $hash{$var}{alleles}[1]) or
                ($ref eq $hash{$var}{alleles}[1] and $alt eq $hash{$var}{alleles}[0])) {
                $hash{$var}{score} = $score || "NA";
            }
        }
    }
    print "[Info] Eigen scores have been added to variants.\n" if $verbose;
    return \%hash
}

##
## Get CADD scores
##
sub get_CADD_GERP {
    my %hash = %{$_[0]};
    my $score = $_[1];

    # Looping throuh all variants and return the Eigen score for all:
    foreach my $var (keys %hash){
        (my $chr = $hash{$var}{GRCh37}[0] ) =~ s/chr//i;
        my $tabix_query = sprintf("tabix %s %s:%s-%s | cut -f1-5,25-28,115-", $caddPath, $chr, $hash{$var}{GRCh37}[2], $hash{$var}{GRCh37}[2]);
        my $lines = `bash -O extglob -c \'$tabix_query\'`;
        #$hash{$var}{CADD} = ["NA", "NA"]; # Initialize CADD scores.
        #$hash{$var}{GERP} = ["NA", "NA", "NA", "NA"]; # Initalize GERP scores
        $hash{$var}{score} = "NA";

        foreach my $line (split("\n", $lines)){
            # Line: 12	56482614	C	C	G	5.52	4.63	510	1.76195e-61	4.931423	25.0
            chomp $line;
            my ($Chrom, $Pos, $ref, $anc, $alt, $GerpN, $GerpS, $GerpRS, $GerpRSpval, $RawScore, $PHRED) = split("\t", $line);

            # Testing alleles:
            if (($ref eq $hash{$var}{alleles}[0] and $alt eq $hash{$var}{alleles}[1]) or
                ($ref eq $hash{$var}{alleles}[1] and $alt eq $hash{$var}{alleles}[0])) {

                    #$hash{$var}{CADD} = [$RawScore, $PHRED]; # Initialize CADD scores.
                    #$hash{$var}{GERP} = [$GerpN, $GerpS, $GerpRS, $GerpRSpval]; # Initalize GERP scores
                    $hash{$var}{score} = $PHRED if $score eq "CADD";
                    $hash{$var}{score} = $GerpS if $score eq "GERP";
            }
        }
    }
    print "[Info] $score scores have been added to variants.\n" if $verbose;
    return \%hash
}

##
## Lifting over the coordinates of the variants to the older build:
##
sub liftover {
    my $gene_name = $_[1];
    my %hash = %{$_[0]};

    my $tempFileName = sprintf($temporaryBedFile, $gene_name);
    open( my $tempbed, ">", $tempFileName) or die "[Error] Temporary bedfile could not be opened for writing: $tempFileName\n";
    foreach my $variant (keys %hash){
        printf $tempbed "%s\t%s\t%s\t%s\n", $hash{$variant}{GRCh38}[0], $hash{$variant}{GRCh38}[1] - 1, $hash{$variant}{GRCh38}[1], $variant;
    }

    # Liftover query:
    my $liftover_query = sprintf("%s %s %s %s_GRCh37.bed %s_unmapped.bed", $liftoverPath, $tempFileName, $chainFile, $gene_name, $gene_name);

    # Add more info if verbose mode is on:
    $liftover_query .= " 2> /dev/null" unless $verbose;

    # Calling liftover:
    `bash -O extglob -c \'$liftover_query\'`;

    # Reading mapped file:
    my $lifted_file = sprintf("%s_GRCh37.bed", $gene_name);
    open(my $lifted, "<", $lifted_file) or die "[Error] After liftover run, the mapped file could not be opened.\n";
    my $liftedVarNo  = 0;
    while (my $line = <$lifted>) {
        chomp $line;
        my ($chr, $start, $end, $SNPID) = split("\t", $line);
        $hash{$SNPID}{"GRCh37"} = [$chr, $start, $end];
        $liftedVarNo ++;
    }

    print "[Info] Number of variants successfully lifted over: $liftedVarNo\n" if $verbose;
    return \%hash;
}


##
## Filtering variation list based on MAF and MAC
##
sub processVar {
    my $variants = $_[0]; # List of all overlapping variants
    my %genotypeContainer = %{$_[1]}; # hash to contain genotype information
    my $total = 0;
    my %hash = ();

    foreach my $variant (split("\n", $variants)){
        $total ++;

        #line: CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	EGAN00001033155
        my ($chr, $pos, $id, $a1, $a2, $qual, $filter, $info, $format, @genotypes) = split(/\t/, $variant);

        # Parsing info field for relevant information:
        $info =~ /AC=(\d+?);.+AN=(\d+?);/;
        my ($ac, $an) = ($1, $2);

        next if length($a2) > 1 or length($a1) > 1; # We don't consider multialleleic sites this time.
        my $MAF = $ac/$an;
        $MAF = 1 - $MAF if $MAF > 0.5;
        my $MAC = $ac;
        $MAC = $an - $ac if $MAC > $an / 2;
        if ($ac >= $Variant{'MAC'} && $MAF <= $Variant{'MAF'}){
            my $SNPID = sprintf("%s:%s_%s/%s", $chr, $pos, $a1, $a2);

            #Storing variant data for all variant:
            $hash{$SNPID}{"alleles"} = [$a1, $a2];
            $hash{$SNPID}{"GRCh38"} = [$chr, $pos];
            $hash{$SNPID}{"frequencies"} = [$ac, $an, $MAF];

            # Parsing genotypes:
            foreach my $gt (@genotypes){
                my @fields = split(":", $gt);
                if    ($fields[0] eq "0/0") { push(@{$genotypeContainer{$SNPID}}, 0) }
                elsif ($fields[0] eq "1/1") { push(@{$genotypeContainer{$SNPID}}, 2) }
                elsif ($fields[0] eq "1/0" or $fields[0] eq "0/1") { push(@{$genotypeContainer{$SNPID}}, 1) }
                else {push(@{$genotypeContainer{$SNPID}}, "-9")}
            }
        }
    }

    print STDERR "[Info] Total number of overlapping variants: $total\n" if $verbose;
    printf STDERR ( "[Info] Number of filtered variants: %s\n", scalar(keys %hash));
    return (\%hash, \%genotypeContainer);
}


##
## bcftools query:
##
sub GetVariants {
    my $merged = $_[0];
    my $variations = '';
    my $distance = 0;

    #
    foreach my $line (split("\n", $merged)){
        my ($chr, $start, $end) = split("\t", $line);
        $distance += $end - $start;
        my $bcftoos_query = sprintf("bcftools query -f \"%%LINE\" %s -r chr%s:%s-%s", $vcfFile, $chr, $start, $end);
        $variations .= `bash -O extglob -c \'$bcftoos_query\'`;
    }

    print sprintf("[Info] Total covered genomic regions: %s bp\n", $distance) if $verbose;
    return $variations;

}

##
## Based on the resulting lines, a new bedfile will be generated
##
sub CollapseRegions {
    my @filteredLines = @_;
    my $bedString = '';

    # Saving files into a bedfile, extend regions if necessary:
    foreach my $line (@filteredLines){
        # chr=3;start=52999670;end=52999869;source=GTEx;
        $line =~ /chr=(.+?);start=(.+?);end=(.+?);.*source=(.+?);.*class=(.+?);/;
        my ($chr, $start, $end, $source, $class) = ($1,$2, $3, $4, $5);
        if ($source eq "GENCODE") {
            $bedString .= sprintf("%s\t%s\t%s\t%s\t%s\n", $chr, $start - $GENCODE{"extend"}, $end + $GENCODE{"extend"}, $source, $class);
        }
        else {
            $bedString .= sprintf("%s\t%s\t%s\t%s\t%s\n", $chr, $start, $end , $source, $class);
        }
    }

    my $queryString = sprintf("mergeBed -i <(echo -e \"%s\")", $bedString);
    my $merged = `bash -O extglob -c \'$queryString\'`;
    return $merged;
}


##
## Filtering lines based on submitted criteria
##
sub FilterLines {
    my $lines = $_[0];
    my @output_lines = ();
    my %hash = ();

    foreach my $line (split (/\n/, $lines)){
        $line =~ /source=(.+?);/;

        # Implementing possible sources:
        next unless $1;
        my $source = $1;

        $line =~ /class=(.+?);/;
        next unless $1;
        my $class = $1;

        if ($source eq "GENCODE" and exists $GENCODE{$class}) {
            next if exists $GENCODE{"minor"} && $line =~ /appris=Minor/;
            push (@output_lines, $line);
            $hash{GENCODE}{$class} ++;
        }
        elsif (($source eq "GTEx" and exists $GTEx{$class})
               or ($source eq "GTEx" and exists $GTEx{'allreg'})){
            push (@output_lines, $line);
            $hash{GTEx}{$class} ++;
        }
        elsif (($source eq "overlap" and exists $overlap{$class})
               or ($source eq "overlap" and exists $overlap{'allreg'})){
            print "$source, $class\n$line\n";
            push (@output_lines, $line);
            $hash{overlap}{$class} ++;

        }
        # It is easy to expand this design.
    }

    if ($verbose) {
        if (exists $hash{GENCODE}) {
            my $report = '';
            foreach my $feature (keys %{$hash{GENCODE}}){
                $report .= "$feature: $hash{GENCODE}{$feature}, ";
            }
            print STDERR "[INFO] From the GENCODE dataset the following features were extracted: $report\n" if $verbose;
        }
        if (exists $hash{GTEx}) {
            my $report = '';
            foreach my $feature (keys %{$hash{GTEx}}){
                $report .= "$feature: $hash{GTEx}{$feature}, ";
            }
            print STDERR "[INFO] The following regulatory features were linked to the gene: $report\n" if $verbose;
        }
        if (exists $hash{overlap}) {
            my $report = '';
            foreach my $feature (keys %{$hash{overlap}}){
                $report .= "$feature: $hash{overlap}{$feature}, ";
            }
            print STDERR "[INFO] The following overlapping regulatory features were extracted: $report\n"  if $verbose;
        }
    }

    print "[Info] Selected lines:\n", join("\n", @output_lines),"\n" if $verbose;

    return @output_lines;
}

##
## This function returns all lines corresponding to the given gene using bedtools.
##
sub BedToolsQuery {
    my ($chr, $start, $end, $stable_ID) = @_;
    my $queryString = sprintf("intersectBed -wb -a <(echo -e \"%s\\t%s\\t%s\\t%s\") -b %s -sorted | grep -w gene_ID=%s | cut -f9-",
                                $chr, $start, $end, $stable_ID, $geneBedFile, $stable_ID);
    print "[Info] IntersectBed query string: $queryString\n" if $verbose;
    my $query = `bash -O extglob -c \'$queryString\'`;
    return $query;

}

##
## This function returns genomic coordinates of a gene given its stable ID or external name.
##
sub GetCoordinates {
    my $ID = $_[0];

    # If Stable ID is given:
    my ($chr, $start, $end, $stable_ID, $name) = ('NA' x 5);
    if ($ID =~ /ENSG/) {
        my $g = $ga->fetch_by_stable_id($ID);
        $stable_ID = $g->stable_id;
        $start = $g->seq_region_start;
        $end = $g->seq_region_end;
        $chr = $g->seq_region_name;
        $name = $g->external_name;
    }
    else { # Assuming gene name is given:
        foreach my $g (@{$ga->fetch_all_by_external_name($ID)}){
            next unless $g->stable_id;
            $stable_ID = $g->stable_id;
            $start = $g->seq_region_start;
            $end = $g->seq_region_end;
            $chr = $g->seq_region_name;
            $name = $g->external_name;
        }
    }

    return ($chr, $start, $end, $stable_ID, $name)
}


##
## Below is the documentation for the script:
##
=pod

=head1 DESCRIPTION

This script was written to select variants for collapsing tests (MONSTER). The script allows to specify
a custom set of criteria to define genomic regions of interest, scoring and weighting methods.

In more details: the selection is perfomed in a gene basis: the user can specify which functional
part of the gene is of interest (eg. exon, CDS stc). Besides the GENCODE elements, Ensembl regulatory
features are also available to be selected: if a feature overlaps with the queried gene, or
overlaps with with an eQTL signal which is linked to the gene.

Overlapping variants from the 15X helic data are then selected using intersectbed. After applying
a various set of filters, scores are assigned to the variants to use them as weights by the
collapsing test. Optionally weigh can be adjusted for allele frequencies to add
more weight to rare variants.

=head1 VERSION

v.1.0. Last modified: 2016.04.11 by Daniel Suveges, mail: ds26@sanger.ac.uk

=head1 SYNOPSIS

perl $0 --input|i <Input file>
        --output|o <Output file prefix>
        --GENCODE|G <feature list>
        --GTEx|E <feature list>
        --overlap|L <featire list>
        --MAF <upper MAF threshold>
        --MAC <lower MAC threshold>
        --extend <bp to extend GENCODE regions in bp>
        --score|s <Scoring method>
        --SkipMinor
        --MAF_weight|w
        --ExpConst|k
        --verbose|v
        --help|h

=head1 OPTIONS AND ARGUMENTS

B<Input file:>

=over 4

List with gene names or Ensembl stable gene IDs (Each ID in new line, separated by a Unix newline character) Required parameter!

=back

B<Output files:>

=over 4

Output prefix ${output}_genotype and ${output}_variants files are saved. Genotype and SNP files are suitable input for MONSTER (more information of the format: I<www.stat.uchicago.edu/~mcpeek/software/MONSTER/MONSTER_v1.2_doc.pdf>) Output prefix is a required parameter!


=back

B<Other parameters>:

=over 4

B<GENCODE>: accepted gencode parameters: comma separated list of GENCODE feature names: gene, exon, transcript, CDS, UTR (without spaces!)

B<GTEx>: specifying the types of Ensembl regulatory features that were linked to the gene by overlapping eQTL signal.

B<Overlap>: specifying the types of Ensembl regulatory features that were overlapping the queried gene.

For GTEx, Overlap the following regulatory feature names are accepted: promoter, CTCF, enhancer, promoterFlank, openChrom, TF_bind, allreg (where allreg equal to selecting all regulatory features.) (without spaces!)

B<MAF>/B<MAC> for setting the minimal allele count and maximal allele frequencies. Both MAF/MAC values are calculated from the corresponding vcf line. (default MAF: 0.05, MAC: 0)

B<score>: available scoring meghods are: GERP, CADD, Eigen (default: Eigen)

B<extend>: the selected GENCODE regions will be extended with the submitted value (default: 0bp)

B<SkipMinor>: any GENCODE feature will be excluded that belongs to a transcript annotated as minor in the APPRIS database.

B<MAF_weight>: MAF weighting will be applied. Exponential scoring the lower MAF regions. C<s_w = s * e^(-k * ($AC - 1)/$AN)>
    where s_w: weighted score, s: score, AC: allele count, AN: allele number, k: exponential constant

B<ExpConst>: describes the steepness of the weighting function. (default: 50)

B<Verbose>: prints out more reports during the progression of the script.

B<help>: prints out this message

=back


=head1 FILES:

To properly run the script all of these files have to be found. Otherwise the script exits:

=over 4

B<$geneBedFile>: GENCODE annotation bedfile. Mapped to GRCh38 coordinates.

B<$vcfFile>: 15x HELIC MANOLIS datafile.

B<$liftoverPath>: As scores are not awailable to the newer build, after getting the variants, we map them to the older build (GRCh37).

B<$chainFile>: chain file to map GRCh38 coordinates to GRCh37 build.

B<$EigenPath>: the genome wide set of eigen scores. It also contains the GERP scores as well.

B<$caddPath>: The genome wide set of CADD scores.

=back
=cut
