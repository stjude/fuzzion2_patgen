package BAMCoverage;
# generate coverage depth, including deletions/soft clips, for
# either a single base or an interval

use strict;
use Exporter;

use Bio::DB::Sam;

use Configurable;
use MiscUtils qw(get_hash_option dump_die);
use ReferenceNameMapper;
use UniqueAlignmentID;

use constant CIGAR_MATCH => "M";
use constant CIGAR_INSERTION => "I";
use constant CIGAR_DELETION => "D";
use constant CIGAR_SOFT_CLIP => "S";
use constant CIGAR_SKIP => "N";
use constant CIGAR_HARD_CLIP => "H";
# see SAM specification.  Maybe these are in Bio modules anywhere?

#my $SAM_QUERY_BUFFER = 15;
my $SAM_QUERY_BUFFER = 50;
# use a fairly large flanking buffer in case we want to use soft-clipped
# regions

@BAMCoverage::ISA = qw(Configurable Exporter);
@BAMCoverage::EXPORT_OK = qw();

use MethodMaker qw(
min_quality
include_deleted_bases
include_soft_clips
min_mapq
allow_optical_pcr_duplicates

bam
rnm
bam_ref
coverage
verbose
		  );

sub new {
  my ($type, %options) = @_;
  my $self = {};
  bless $self, $type;
  $self->min_quality(15);
  # as used by standard coverage code
  $self->min_mapq(0);
  # required for RNA, also used by basic Ace2.SAMFinneyCoverage
  $self->allow_optical_pcr_duplicates(0);
  # used by basic Ace2.SAMFinneyCoverage
  # HOWEVER this may not be desirable for some RNA projects!
  $self->include_deleted_bases(0);
  # basic Ace2.SAMFinneyCoverage uses aligned blocks only, so
  # deleted regions are NOT counted.
  # pro: deletions detectable by coverage
  # con: relative ratio of reads w/deletion to total coverage off vs. SNVs?

#  $self->verbose(1);
  $self->configure(%options);
  $self->setup();
  return $self;
}

sub setup {
  my ($self, %options) = @_;

  my $f_bam = $self->bam || die "-bam";

  my $bam = Bio::DB::Sam->new(
			      "-bam" => $f_bam,
			      "-expand_flags" => 1,
			     );
  my $rnm = new ReferenceNameMapper();
  foreach my $chrom ($bam->seq_ids()) {
    # standardize reference names
    $rnm->add_name($chrom);
  }
  $self->rnm($rnm);
  $self->bam_ref($bam);
}

sub get_coverage {
  my ($self, %options) = @_;
  my $chr_raw = get_hash_option(\%options, "-chr");
  my $q_pos = get_hash_option(\%options, "-pos");
  my $q_end = $options{"-end"};
  # optional end position for coverage of 

  my $q_chr = $self->rnm->find_name($chr_raw) || die "no BAM chrom for $chr_raw";
#  printf STDERR "start coverage for %s %s %d\n", $self->bam, $chr_raw, $q_pos;

  my $qs = $q_pos - $SAM_QUERY_BUFFER;
  my $qe = ($q_end || $q_pos) + $SAM_QUERY_BUFFER;
  my $verbose = $self->verbose();

  my @alignments = $self->bam_ref->get_features_by_location(
						  -seq_id => $q_chr,
						  -start  => $qs,
						  -end    => $qe
						 );
  my %total_coverage;
  my %coverage2id;

  $self->coverage(\%total_coverage);
  my $min_mapq = $self->min_mapq();
  my $allow_optical_pcr_duplicates = $self->allow_optical_pcr_duplicates();
  my $min_quality = $self->min_quality();
  my $include_del = $self->include_deleted_bases();
  my $include_soft_clip = $self->include_soft_clips();
  my $uid = new UniqueAlignmentID();

  if ($verbose) {
    printf STDERR "raw count:%d\n", scalar @alignments;
    printf STDERR "include softclip: %d\n", $include_soft_clip || 0;
  }

  foreach my $a (@alignments) {
    my $cigar_a = $a->cigar_array();
    my $usable = 1;
    $usable = 0 if $a->qual() < $min_mapq;
    # minimum mapping quality

    unless ($allow_optical_pcr_duplicates) {
      my $optical_pcr_dup = $a->get_tag_values("DUPLICATE");
      die "optical/pcr dup not defined" unless defined $optical_pcr_dup;
      $usable = 0 if $optical_pcr_dup;
    }

    if ($usable) {
      #
      #  Walk through CIGAR and track.
      #
      my $ref_pos = $a->start;
      # reference base alignment start (1-based)
      my $query_idx = $a->query->start - 1;
      # query sequence index (0-based)

      my $unique_read_id = $uid->get_unique_id($a) if $q_end;

      my $query_dna = $a->query->dna();
      my $query_qual = $a->qscore();
      die "ERROR: seq/quality length mismatch" unless length($query_dna) == @{$query_qual};

      printf STDERR "scanning %s strand:%s raw start=%d seq=%s cigar=%s qi=%d\n", $a->query->seq_id(), $a->strand, $ref_pos, $query_dna, $a->cigar_str, $query_idx if $verbose;

      for (my $i = 0; $i < @{$cigar_a}; $i++) {
	my $c = $cigar_a->[$i];
	# walk through CIGAR
	my ($ctype, $clen) = @{$c};
	printf STDERR "CIGAR idx %d: %s %s\n", $i, $ctype, $clen if $verbose;

	if ($ctype eq CIGAR_SOFT_CLIP) {
	  printf STDERR "start soft $include_soft_clip\n" if $verbose;
	  if ($include_soft_clip) {
	    if ($i == 0) {
	      # leading soft-clips are ignored in query start
	      $ref_pos -= $clen;
	      $query_idx -= $clen;
	      printf STDERR "adjusting ref start to %d\n", $ref_pos if $verbose;
	    }
	    printf STDERR "in soft\n" if $verbose;
	    my $q_qual;
	    for (my $j = 0; $j < $clen; $j++) {
	      $q_qual = $query_qual->[$query_idx];
	      if ($q_qual >= $min_quality) {
		$coverage2id{$ref_pos}{$unique_read_id} = 1 if $q_end;
		$total_coverage{$ref_pos}++;
		printf STDERR "covg %d softclip for %s %s\n", $ref_pos, $self->bam, $a->query->seq_id() if $verbose;
	      }
	      $ref_pos++;
	      $query_idx++;
	    }
	  }
	} elsif ($ctype eq CIGAR_MATCH) {
	  # match or mismatch: affects both query and reference
	  my $q_qual;
	  for (my $j = 0; $j < $clen; $j++) {
	    $q_qual = $query_qual->[$query_idx];
	    if ($q_qual >= $min_quality) {
	      $coverage2id{$ref_pos}{$unique_read_id} = 1 if $q_end;
	      $total_coverage{$ref_pos}++;
	      printf STDERR "covg %d main for %s %s\n", $ref_pos, $self->bam, $a->query->seq_id() if $verbose;
	    } else {
	      printf STDERR "covg %d main_quality_fail %d for %s %s\n", $ref_pos, $q_qual, $self->bam, $a->query->seq_id() if $verbose;
	    }

	    $ref_pos++;
	    $query_idx++;
	  }
	} elsif ($ctype eq CIGAR_SKIP) {
	  # skip: affects the reference but not the query
	  $ref_pos += $clen;
	} elsif ($ctype eq CIGAR_INSERTION) {
	  # insertion: affects query but not reference
	  $query_idx += $clen;
	} elsif ($ctype eq CIGAR_DELETION) {
	  # deletion: affects reference but not query
	  for (my $j = 0; $j < $clen; $j++) {
	    if ($include_del) {
	      $coverage2id{$ref_pos}{$unique_read_id} = 1 if $q_end;
	      $total_coverage{$ref_pos}++;
	      printf STDERR "covg %d del for %s %s\n", $ref_pos, $self->bam, $a->query->seq_id() if $verbose;
	    }

	    $ref_pos++;
	  }
	} elsif ($ctype eq CIGAR_HARD_CLIP) {
	  # hard clip: clipped sequence is NOT present in sequence.
	  #   - leading clips are NOT reflected in alignment start,
	  #     i.e. position is already past the hard clip.
	  #   - trailing clips are ignored anyway.
	  # in essence we can ignore
	  printf STDERR "covg %d hardclip for %s %s (ignored)\n", $ref_pos, $self->bam, $a->query->seq_id() if $verbose;
	  $query_idx -= $clen if $i == 0;
	  # leading hard clip: query pointer incorrect in these cases,
	  # seems to imply hard-clipped sequence is present in read, it isn't
	} else {
	  die sprintf "unhandled CIGAR entry %s read=%s start=%d", $ctype, $a->query->seq_id(), $a->start;
	}
      }
    }
  }

  if ($q_end) {
    # multi-base coverage, e.g. deletion/complex
    my %all;
    for (my $i = $q_pos; $i <= $q_end; $i++) {
      foreach my $id (keys %{$coverage2id{$i}}) {
	unless ($all{$id}) {
#	  printf STDERR "hit %s at %d\n", $id, $i;
	  $all{$id} = 1;
	}
      }
    }

    return scalar keys %all;
  } else {
    # single-base coverage
    return $total_coverage{$q_pos} || 0;
    # can be zero if no reads
  }

}

1;

# LRF support:
# _______________      ______________________________      _______________
#                \____/                              \____/
