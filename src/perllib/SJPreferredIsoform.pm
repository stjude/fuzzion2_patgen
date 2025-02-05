package SJPreferredIsoform;
# parse/retrieve SJ preferred isoforms, e.g.
# /nfs_exports/apps/gnu-apps/NextGen/nextgensupport2/NCBI/RefSeq/gene_transcript_matrix.withNM.mod
# TO DO:
# - look up file from genome config if not specified

use strict;
use Carp qw(confess);

use Configurable;
use Exporter;
use WorkingFile;
use GeneSymbolMapper qw(new_gsm_lite);

@SJPreferredIsoform::ISA = qw(Configurable Exporter);
@SJPreferredIsoform::EXPORT_OK = qw();

use MethodMaker qw(
	file
gene2preferred
nm2preferred
nm2rank
nm2gene
gene2nms
genes_ordered
gene2lines

auto_gsm
gsm
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
  my ($self) = @_;
  my $file = $self->file || die "-file";
  open(SJPI, $file) || die "can't open $file: $!";
  my %preferred;
  my %by_nm;
  my %gene2nms;
  my @genes_ordered;
  my %gene2lines;
  my %nm2rank;
  my %nm2gene;
  my %rank;
  if ($self->auto_gsm()) {
    $self->gsm(new_gsm_lite());
    # 1/2020: not sure why this is on by default and required an
    # object pass-in, possibly in case older code casually embedded
    # SJPI instntiations in loops or something.
  } elsif (!$self->gsm()) {
    # if not already manually provided by user
    printf STDERR "WARNING: SJPI auto gene symbol mapping NOT enabled, consider auto_gsm in constructor\n";
    # since this is a reasonable "do what I mean option" seems
    # OK to introduce a nag
  }

  my $gsm = $self->gsm();
  while (my $line = <SJPI>) {
    chomp $line;
    $line =~ s/\r$//;
    my @f = split /\t/, $line;
    die unless @f == 3;
    my ($gene, $db_id, $nm) = @f;
    if ($gsm) {
      $gsm->add_gene("-gene" => $gene);
    }

    my $primary;
    $nm2gene{$nm} = $gene;
    if ($preferred{$gene}) {
      $primary = 0;
    } else {
      $primary = 1;
      $preferred{$gene} = $nm;
      push @genes_ordered, $gene;
    }
    my $rank = $rank{$gene} || 1;
    $nm2rank{$nm} = $rank;
    $rank{$gene} = $rank + 1;

#    printf STDERR "by_nm %s %s %s\n", $gene, $nm, $primary;

    $by_nm{$nm} = $primary;
    push @{$gene2nms{$gene}}, $nm;
    push @{$gene2lines{$gene}}, $line;
  }
  close SJPI;
  $self->gene2preferred(\%preferred);
  $self->nm2preferred(\%by_nm);
  $self->nm2rank(\%nm2rank);
  $self->nm2gene(\%nm2gene);
  $self->gene2nms(\%gene2nms);
  $self->gene2lines(\%gene2lines);
  $self->genes_ordered(\@genes_ordered);
}

sub get_preferred_isoform {
  my ($self, $gene) = @_;
  my $g2p = $self->gene2preferred();
  my $idx;
  if ($g2p->{$gene}) {
    $idx = $gene;
  } elsif (my $gsm = $self->gsm) {
    # attempt to resolve gene symbols, e.g. for KMT2C the SJPI
    # entry is for old symbol MLL3
    $idx = $gsm->resolve("-symbol" => $gene) || $gene;
  }

  return $g2p->{$idx};
}

sub is_preferred_nm {
  my ($self, $nm) = @_;
  $nm = strip_version($nm);
  return $self->nm2preferred()->{$nm};
}

sub get_nm_count {
  my ($self, $gene) = @_;
  return scalar @{$self->gene2nms->{$gene}};
}

sub get_accessions {
  my ($self, $gene) = @_;
  return $self->gene2nms->{$gene};
}

sub write_preferred_file {
  my ($self, %options) = @_;
  my $reorder = $options{"-reorder"} || {};
  my $outfile = "sjpi.tab";

  my $wf = new WorkingFile($outfile);
  my $fh = $wf->output_filehandle();
  foreach my $gene (@{$self->genes_ordered}) {
    my $lines = $self->gene2lines->{$gene} || die;
    my @postponed;
    my $promote = $reorder->{$gene};
    my @write;
    if ($promote) {
      my ($first, @rest);
      foreach my $l (@{$lines}) {
	my @f = split /\t/, $l;
	my $acc = $f[2];
	if ($acc eq $promote) {
	  $first = $l;
	} else {
	  push @rest, $l;
	}
      }
      die unless $first and @rest;
      @write = ($first, @rest);
#      die join "\n", @{$lines}, "---", @write;
    } else {
      @write = @{$lines};
    }

    foreach my $line (@write) {
      printf $fh "%s\n", $line;
    }

  }
  $wf->finish;
}

sub get_genes {
  my ($self) = @_;
  return [ sort keys %{$self->gene2nms} ];
}

sub override_preferred_isoform {
  my ($self, %options) = @_;
  my $gene = $options{"-gene"} || die;
  my $wanted_nm = $options{"-accession"} || die;
  my $all_nm = $self->gene2nms()->{$gene} || die "no accessions for $gene";
  die "can't find $wanted_nm in list for $gene" unless grep {$_ eq $wanted_nm} @{$all_nm};

  foreach my $nm (@{$all_nm}) {
    my $is_preferred = $nm eq $wanted_nm ? 1 : 0;
    $self->nm2preferred->{$nm} = $is_preferred;
  }
}

sub get_nm_rank {
  my ($self, $nm, $gene) = @_;
  $nm = strip_version($nm);
  my $rank = $self->nm2rank()->{$nm};
  if ($gene and my $nm_list = $self->gene2nms->{$gene}) {
    # optional feature: use user-provided gene symbol to indicate preference.
    # e.g. 1.147575870.A.G, NBPF11, VEP annotates NM_001101663
    # and NM_183372 but both of these are ranked 1 in SJPI
    # because SJPI puts NM_183372 in a separate list for NBPF24,
    # which is a previous symbol for NBPF11.
    $rank = undef unless grep {$_ eq $nm} @{$nm_list};
  }
  return $rank;
}

sub get_nm_gene {
  my ($self, $nm) = @_;
  return $self->nm2gene()->{$nm};
}

sub get_gene2nms {
  my ($self) = @_;
  return $self->gene2nms();
}

sub get_genes_ordered {
  my ($self) = @_;
  return $self->genes_ordered();
}

sub get_lines {
  my ($self, $gene) = @_;
  return $self->gene2lines->{$gene};
}

sub strip_version {
  # STATIC
  my ($acc) = @_;
  $acc =~ s/\.\d+$//;
  return $acc;
}

1;

# LRF support:
# _______________      ______________________________      _______________
#                \____/                              \____/               
