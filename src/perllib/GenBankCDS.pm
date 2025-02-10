package GenBankCDS;
# GenBank coding sequence tools

use strict;
use Exporter;
use Carp qw(confess);

use Bio::SeqIO;
use Bio::Tools::CodonTable;

use Configurable;
use MiscUtils qw(get_hash_option dump_die);
use FileUtils qw(universal_open);

@GenBankCDS::ISA = qw(Configurable Exporter);
@GenBankCDS::EXPORT_OK = qw(tolerant_protein_compare);

use MethodMaker qw(
		    file
		    genbank_info
		    die_if_not_found
		    is_ambiguous

cds_start
cds_end
base2exon_number
		  );

sub new {
  my ($type, %options) = @_;
  my $self = {};
  bless $self, $type;
  $self->configure(%options);
  $self->setup(%options);
  return $self;
}

sub setup {
  my ($self, %options) = @_;
  my $f_gb = $self->file() || die "-file";
  # file containing a single GenBank record
  # TO DO: accept record object directly?

  my $fh = universal_open($f_gb);
  my $stream = Bio::SeqIO->new("-fh" => $fh,
			       -format => 'GenBank');
  my $record = $stream->next_seq();
  # Bio::SeqIO::genbank
  my $info = $self->parse_one_genbank($record);

  $self->genbank_info($info);
}

sub parse_one_genbank {
  # STATIC
  # Bio::SeqIO::genbank
  # derived from fz2/fusion2breakpoints.pl
  my ($self, $record) = @_;
  my $verbose = 0;

  my $accession = $record->accession_number();
  my $version = $record->version();

  my $accession_versioned = join ".", $accession, $version;
  # keep version
  my $desc = $record->desc() || die "no desc for $accession";

  my $cds_count = 0;

  my %r;
  $r{accession} = $accession;
  $r{version} = $version;
  $r{accession_versioned} = $accession_versioned;

  my %base2exon_number;
  # NOTE: refers to "mrna_full", *NOT* CDS

  print STDERR "\n" if $verbose > 1;
  printf STDERR "starting %s...\n", $accession if $verbose;

  my $exon_counter = 0;

  my %genes;
  foreach my $feature ($record->get_SeqFeatures()) {
    printf STDERR "  %s\n", $feature->primary_tag() if $verbose > 1;
    my $ptag = $feature->primary_tag();
    if ($ptag eq "gene") {
      my @values = $feature->get_tag_values($ptag);
      die if @values > 1;
      $genes{$values[0]} = 1;
    } elsif ($ptag eq "CDS") {
      $cds_count++;
      foreach my $tag ($feature->get_all_tags()) {
	printf STDERR "tag:%s values:%s\n", $tag, join ",", $feature->get_tag_values($tag) if $verbose > 1;
	if ($tag eq "translation") {
	  my @values = $feature->get_tag_values($tag);
	  die unless @values == 1;
	  my $protein = $values[0];
	  $r{protein} = $protein;

	  my $cds_start = $feature->start() || die;
	  my $cds_end = $feature->end() || die;
	  $self->cds_start($cds_start);
	  $self->cds_end($cds_end);

	  my $mrna = $record->seq() || die;

	  my $utr5_plus_cds = substr($mrna, 0, $cds_end);

	  $cds_end -= 3;
	  # trim off stop codon for purposes of regenerating coding sequence

	  my $mrna_cds = substr($mrna, $cds_start - 1, ($cds_end - $cds_start) + 1);
	  my $protein_check = mrna2protein($mrna_cds);
#	  printf STDERR "cds start:%d end:%d\n", $cds_start, $cds_end;
#	  printf STDERR "mrna len: %d\n", length $mrna;
#	  printf STDERR "mrna cds len: %d\n", length $mrna_cds;

	  die sprintf "ERROR: protein mismatch parsing GenBank record for %s:\n%s\n%s\n", $accession, $protein, $protein_check unless tolerant_protein_compare($protein, $protein_check);
	  # QC check: verify RNA for CDS produces annotated protein

	  $r{mrna_cds} = $mrna_cds;
	  $r{mrna_cds_with_utr5} = $utr5_plus_cds;
	  $r{mrna_full} = $mrna;
	  # unedited, required for using base2exon_number!
	}
      }
    } elsif ($ptag eq "exon") {
      my $exon_number = ++$exon_counter;
      my $start = $feature->start;
      my $end = $feature->end;
      for (my $pos = $start; $pos <= $end; $pos++) {
	$base2exon_number{$pos} = $exon_number;
      }
    }
  }
  $r{gene} = join ",", sort keys %genes;

  die "no info!" unless %r;
#  die sprintf "CDS count %d for %s", $cds_count, $accession if $cds_count != 1;
  die "multiple CDS" if $cds_count > 1;
  # will need to return multiple records if this happens

  $self->base2exon_number(\%base2exon_number);

  return \%r;
}

sub mrna2protein {
  # STATIC
  # TO DO: move to misc utilities code?
  my ($mrna) = @_;
  my $ct = new Bio::Tools::CodonTable();
  my @codons = $mrna =~ /(...)/g;
  my $protein = join "", map {$ct->translate($_)} @codons;
  return $protein;
}

sub tolerant_protein_compare {
  # STATIC
  my ($p1, $p2) = @_;
  my $ok;

  if (length($p1) == length($p2)) {
    $ok = 1 if $p1 eq $p2;
    # perfect match

    unless ($ok) {
      # check for alternative initiation codon, e.g.:
      # https://www.ncbi.nlm.nih.gov/nuccore/NM_001354870.1/
      #
      # CDS             1161..2522
      #         /gene="MYC"
      #         /gene_synonym="bHLHe39; c-Myc; MRTL; MYCC"
      #         /note="isoform 2 is encoded by transcript variant 2;
      #         non-AUG (CUG) translation initiation codon;
      if (substr($p1,0,1) eq "M" or substr($p2,0,1) eq "M") {
	# one of the sequences being compared starts with a canonical
	# M (in is example, refSeq sequence shows start codon of M
	# even though "CTG" at position 1161 encodes L instad:
#
#	1141 tccttgcagc tgcttagacg ctggattttt ttcgggtagt ggaaaaccag cctcccgcga
#
	my $p1_alt = substr($p1, 1);
	my $p2_alt = substr($p2, 1);
	$ok = 1 if $p1_alt eq $p2_alt;
      }
    }

    unless ($ok) {
      # exceptions to stop codon translation:
      #
      # https://www.ncbi.nlm.nih.gov/nuccore/NM_001315537.2
      #   stop codon readthrough signal

      # https://www.ncbi.nlm.nih.gov/nuccore/NM_182743.3:
      #
      # "The 3' UTRs of selenoprotein mRNAs contain a conserved
      # stem-loop structure, the Sec insertion sequence (SECIS)
      # element, which is necessary for the recognition of UGA as a
      # Sec codon rather than as a stop signal."
      my $len = length($p1);
      my @mismatch_i;
      for (my $i = 0; $i < $len; $i++) {
	push @mismatch_i, $i unless substr($p1, $i, 1) eq substr($p2, $i, 1);
      }
      if (@mismatch_i == 1) {
	my $i = $mismatch_i[0];
	my %both = map {$_, 1} (substr($p1, $i, 1),
				substr($p2, $i, 1)
			       );
	if ($both{"*"} and ($both{"X"} or $both{"U"})) {
	  # NM_001315537, GenBank sequence shows X instead of *
	  # NM_182743.3: GenBank shows U instead of *
	  $ok = 1;
	}
      }
    }
  } else {
    $ok = 0;
  }
  return $ok;
}

sub get_rna_for_protein {
  # get equivalent RNA for a protein subsequence
  my ($self, $pchunk) = @_;
  my $gbi = $self->genbank_info() || die;

  my $protein = $gbi->{protein};
  my $rna_full = $gbi->{mrna_cds} || die;

  my @hits = $protein =~ /$pchunk/ig;
  my $rna = "";
  my $is_ambiguous = 0;
  if (@hits == 1) {
    # OK
    my $idx = index(uc($protein), uc($pchunk));
    die "lookup failed" if $idx == -1;
    $rna = substr($rna_full, $idx * 3, length($pchunk) * 3);
    die "protein sanity check fail" unless mrna2protein($rna) eq $pchunk;
  } elsif (@hits > 1) {
    # given sequence ambiguous in the protein
    $rna = "";
    $is_ambiguous = 1;
  } elsif ($self->die_if_not_found()) {
    # shouldn't be the default, could be e.g. searching for a sequence
    # that isn't present in a particular sub-isoform
    confess sprintf "fail: %d hits for %s in %s\n",
      scalar(@hits), $pchunk, $protein;
  }
  $self->is_ambiguous($is_ambiguous);

  return $rna;
}

sub get_info {
  my ($self) = @_;
  return $self->genbank_info();
}

sub get_feature_for_base_number {
  my ($self, $bn) = @_;

  my $exno = $self->base2exon_number()->{$bn} || die "can't find exon number for base # $bn";
  my $fname;
  if ($bn < $self->cds_start()) {
    $fname = "5utr";
  } elsif ($bn <= $self->cds_end()) {
    $fname = "exon";
  } else {
    $fname = "3utr";
  }
  return sprintf "%s_%d", $fname, $exno;
}

1;

# LRF support:
# _______________      ______________________________      _______________
#                \____/                              \____/
