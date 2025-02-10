package GeneSymbolStandardizer;
# standardize gene symbols based on NCBI gene_info.table
# ftp://ftp.ncbi.nih.gov/gene/DATA/README
# ftp://ftp.ncbi.nih.gov/gene/DATA/gene_info.gz
# Michael Edmonson 12/2013

# needs work for ambiguous synonyms: e.g. HNPCC maps to multiple genes!

use strict;
use Exporter;

use Configurable;
use DelimitedFile;

@GeneSymbolStandardizer::ISA = qw(Configurable Exporter);
@GeneSymbolStandardizer::EXPORT_OK = qw();

use MethodMaker qw(
                   filename
                   gene_canonical
                   gene_synonym
		  );

sub new {
  my ($type, %options) = @_;
  my $self = {};
  bless $self, $type;
  $self->configure(%options);
  $self->setup();
  return $self;
}

sub setup {
  my ($self, %options) = @_;

  my $fn = $self->filename || die "-filename";
  my $df = new DelimitedFile(
			  "-file" => $fn,
			  "-delimiter" => "\t",
			  "-headers" => 1,
			 );
  # hack:
  # 1. filtered to tax_id 9606 (human)
  # 2. tabified header line

  my %gene_canonical;
  my %gene_synonym;
  # this should only be used secondarily, e.g. PROC maps to APC [fail]
  while (my $row = $df->get_hash()) {
    my $ncbi_sym = $row->{Symbol};
    my $canonical = $row->{Symbol_from_nomenclature_authority} || die;
    my @other;
    if ($canonical eq "-") {
      $canonical = $ncbi_sym;
    } elsif ($ncbi_sym ne $canonical) {
      # treat as synonym
      push @other, $ncbi_sym;
    }
    $gene_canonical{$canonical} = 1;

    push @other, split /\|/, $row->{Synonyms};
    foreach my $s (@other) {
#      printf STDERR "secondary %s => %s\n", $s, $canonical;
      die if $s eq $canonical;
      $gene_synonym{$s}{$canonical} = 1;
      # might be ambiguous
    }
  }

  foreach my $s (keys %gene_canonical) {
#    printf STDERR "%s found in synonym list (maps to %s)\n", $s, $gene_synonym{$s} if $gene_synonym{$s};
    delete $gene_synonym{$s};
    # delete "synonyms" which are actually canonical symbols.
    # e.g. APC is supposedly a synonym for PROC [fail].
  }

  $self->gene_canonical(\%gene_canonical);
  $self->gene_synonym(\%gene_synonym);
}

sub is_canonical {
  my ($self, $symbol, %options) = @_;
  return $self->gene_canonical->{$symbol};
}

sub find {
  my ($self, $symbol, %options) = @_;
  my $result;
  if ($self->gene_canonical->{$symbol}) {
    $result = $symbol;
  } elsif (my $to = $self->gene_synonym->{$symbol}) {
    $result = [ keys %{$to} ];
  }
  return $result;
}

1;

# LRF support:
# _______________      ______________________________      _______________
#                \____/                              \____/               
