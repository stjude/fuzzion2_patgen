package TranscriptFlanking;
# get upstream/downstream transcript sequence via refFlat mapping
# MNE 3/2021

use strict;
use Exporter;
use Carp qw(confess);

use Configurable;
use MiscUtils qw(get_hash_option dump_die);
use GenomeUtils qw(reverse_complement);

@TranscriptFlanking::ISA = qw(Configurable Exporter);
@TranscriptFlanking::EXPORT_OK = qw();

use MethodMaker qw(
	fai
		  );

sub new {
  my ($type, %options) = @_;
  my $self = {};
  bless $self, $type;
  $self->configure(%options);
  return $self;
}

sub get_upstream_downstream {
  my ($self, %options) = @_;
  # - upstream sequence is from start of the gene model through given base
  # - downstream sequence is from given base to end of gene model
  # - returned sequence is native to the transcript strand
  #
  # TO DO:
  # - policy for UTR sequence?
  my $row = $options{"-row"} || die;
  my $pos = $options{"-pos"} || die;

  confess sprintf 'pos "%s" is not numeric', $pos unless $pos =~ /^\d+$/;
  # in case space characters, etc. creep in from source spreadsheets

  my $dir = $options{"-direction"} || die;
  my $fai = $self->fai() || die "-fai";

  my $pos_only = $options{"-pos-only"};
  # optional

  my $strand = $row->{strand} || die;
  my $chrom = $row->{chrom} || die;
  my $result;

  my $exons = $row->{exons};
  if (1) {
    printf STDERR "raw pos: %s.%d\n", $chrom, $pos;
    printf STDERR "exons: %s\n", join ", ", map {sprintf "%d-%d", @{$_}{qw(start end)}} @{$exons};
  }

  my $ebi = find_exon_block_index($row, $pos);

  unless (defined $ebi) {
    # if the position falls outside an exon, try to find the nearest
    # exon in the appropriate stream direction

    my $new_pos;
    for (my $i = 0; $i < @{$exons} - 1; $i++) {
      if ($pos >= $exons->[$i]->{end} and
	  $pos <= $exons->[$i + 1]->{start}) {
	# position is intronic
	if ($strand eq "+") {
	  if ($dir eq "downstream") {
	    $new_pos = $exons->[$i + 1]->{start};
	    last;
	  } elsif ($dir eq "upstream") {
	    $new_pos = $exons->[$i]->{end};
	  } else {
	    die;
	  }
	} elsif ($strand eq "-") {
	  if ($dir eq "downstream") {
	    $new_pos = $exons->[$i]->{end};
	  } elsif ($dir eq "upstream") {
	    $new_pos = $exons->[$i + 1]->{start};
	  } else {
	    die;
	  }
	} else {
	  die;
	}
      }
    }

    if ($new_pos) {
      $pos = $new_pos;
      printf STDERR "  intron-to-exon pos: %d\n", $pos;
      $ebi = find_exon_block_index($row, $pos);
      die "intron-to-exon map failed" unless defined $ebi;
    }
  }

  return $pos if $pos_only;

#  printf STDERR "  adj exon pos: %s\n", $ebi || "undef";

  if (defined $ebi) {
    # position is within an exon interval
    if ($strand eq "+") {
      if ($dir eq "downstream") {
	$result = genomic_down($exons, $ebi, $chrom, $pos, $fai);
      } elsif ($dir eq "upstream") {
	$result = genomic_up($exons, $ebi, $chrom, $pos, $fai);
      } else {
	die $dir;
      }
    } elsif ($strand eq "-") {
      if ($dir eq "downstream") {
	$result = genomic_up($exons, $ebi, $chrom, $pos, $fai);
      } elsif ($dir eq "upstream") {
	$result = genomic_down($exons, $ebi, $chrom, $pos, $fai);
      } else {
	die;
      }
      $result = reverse_complement($result);
    } else {
      die;
    }
  }

  if ($options{"-extended"}) {
    die unless wantarray;
    return ($result, $pos);
  } else {
    return $result;
  }
}

sub find_exon_block_index {
  # STATIC
  # TO DO: option if intronic to find nearest boundary
  my ($row, $pos) = @_;
  my $ei;
  my $exons = $row->{exons};
  for ($ei = 0; $ei < @{$exons}; $ei++) {
#    printf STDERR "%d-%d\n", @{$exons->[$ei]}{qw(start end)};
    if ($pos >= $exons->[$ei]->{start} and
	$pos <= $exons->[$ei]->{end}) {
      return $ei;
    }
  }
  return undef;
}

sub genomic_down {
  # STATIC
  my ($exons, $ebi, $chrom, $pos, $fai) = @_;

  my @blocks;
  for (my $i = $ebi; $i < @{$exons}; $i++) {
    my ($start, $end) = @{$exons->[$i]}{qw(start end)};
    $start = $pos if $i == $ebi;
    push @blocks, [ $start, $end ];
  }
  printf STDERR "blocks: %s\n", join ", ", map {join "-", @{$_}} @blocks;

  my $result = "";
  foreach my $block (@blocks) {
    my ($start, $end) = @{$block};
    my $chunk = $fai->get_chunk(
				"-id" => $chrom,
				"-start" => $start,
				"-end" => $end
			       );
    printf STDERR "%d-%d: %s\n", $start, $end, $chunk;
    $result .= $chunk;
  }
  return $result;
}

sub genomic_up {
  # STATIC
  my ($exons, $ebi, $chrom, $pos, $fai) = @_;

  my @blocks;
  for (my $i = $ebi; $i >= 0; $i--) {
    my ($start, $end) = @{$exons->[$i]}{qw(start end)};
    $end = $pos if $i == $ebi;
    push @blocks, [ $start, $end ];
  }
  printf STDERR "blocks: %s\n", join ", ", map {join "-", @{$_}} @blocks;

  my $result = "";
  foreach my $block (@blocks) {
    my ($start, $end) = @{$block};
    my $chunk = $fai->get_chunk(
				"-id" => $chrom,
				"-start" => $start,
				"-end" => $end
			       );
    printf STDERR "%d-%d: %s\n", $start, $end, $chunk;
    $result = $chunk . $result;
  }
  return $result;
}

1;

# LRF support:
# _______________      ______________________________      _______________
#                \____/                              \____/
