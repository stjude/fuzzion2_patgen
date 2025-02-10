package ReadthroughFusionAnnotator;
# given a gene pair (and optionally a contig sequence), predict
# whether event is a readthrough fusion
# MNE 2/2022

use strict;
use Exporter;

use Configurable;
use MiscUtils qw(get_hash_option dump_die);
use GeneSymbolMapper qw(new_gsm_lite);
use TemporaryFileWrangler;
use BLATClient;
use File::Copy;
use GenomeUtils qw(cook_chromosome_name);

@ReadthroughFusionAnnotator::ISA = qw(Configurable Exporter);
@ReadthroughFusionAnnotator::EXPORT_OK = qw();

my $READTHROUGH_MAX_DISTANCE_NEAR = 200000;
# canonical / more likely
my $READTHROUGH_MAX_DISTANCE_FAR = 1000000;
# further / less likely

use MethodMaker qw(
        rf
	max_distance_near
	max_distance_far
	gsm
	rt_interval
	rt_info
		  );
# rf = RefFlatFile.pm

sub new {
  my ($type, %options) = @_;
  my $self = {};
  bless $self, $type;
  $self->max_distance_near($READTHROUGH_MAX_DISTANCE_NEAR);
  $self->max_distance_far($READTHROUGH_MAX_DISTANCE_FAR);
  $self->configure(%options);
  $self->setup();
  return $self;
}

sub setup {
  my ($self) = @_;
  my $rf = $self->rf() || die "RefFlatFile not initialized";

  # init refFlat gene symbol disambiguation:
  my $gsm = new_gsm_lite();
  foreach my $r (@{$rf->rows}) {
    my $gene = $r->{gene} || die;
    $gsm->add_gene("-gene" => $gene);
  }
  $self->gsm($gsm);
}

sub is_readthrough_fusion {
  my ($self, %options) = @_;
  my $geneA_raw = $options{"-geneA"} || die "-geneA";
  my $geneB_raw = $options{"-geneB"} || die "-geneB";
  my $contig = $options{"-contig"};
  # to be a potential readthrough candidate, genes must:
  # - have the same chromosome and strand
  # - be consecutive on the strand, with A coming before B
  # - be within a max distance (default 200k nt)
  my $genome = $options{"-genome"};
  my $blat_host = $options{"-blat-host"};
  my $blat_port = $options{"-blat-port"};

  my $bc = new BLATClient();
  $bc->outfile_literal(1);

  my $verbose = $ENV{VERBOSE};

  my @bco;
  if ($contig) {
    # fusion contig sequence to be blat'd
    printf STDERR "contig: %s\n", $contig if $verbose;
    if ($genome) {
      $bc->config_genome($genome);
    } elsif ($blat_host and $blat_port) {
      $bc->host($blat_host);
      $bc->port($blat_port);
    } else {
      die "-contig requires -genome GENOME | -blat-host HOST -blat-port PORT";
    }
  }

  my $rf = $self->rf() || die;
  my $gsm = $self->gsm() || die;
  my $max_distance_near = $self->max_distance_near() || die;
  my $max_distance_far = $self->max_distance_far() || die;

  my @rt_rf_pairs;

  my $has_distance_near;
  my $has_distance_far;
  my $has_consistent_strand;
  my $has_a_b_ordering;
  my $has_second_to_last_to_second;
  my $has_gene_overlap;

  if ($geneA_raw ne $geneB_raw) {
    my $geneA_rf = $gsm->find($geneA_raw);
    my $geneB_rf = $gsm->find($geneB_raw);

    if ($geneA_rf and $geneB_rf) {
      my $rowsA = $rf->find_by_gene($geneA_rf) || die;
      my $rowsB = $rf->find_by_gene($geneB_rf) || die;

      foreach my $ra (@{$rowsA}) {
	my $ra_chrom = $ra->{chrom} || die;
	my $ra_strand = $ra->{strand} || die;

	foreach my $rb (@{$rowsB}) {
	  if ($ra_chrom eq $rb->{chrom}) {
	    # - gene mappings are on same chromosome (there may be be
	    #   multiple chrom mappings)

	    my $pair_strand_ok = $ra_strand eq $rb->{strand};

	    $has_consistent_strand = 1 if $pair_strand_ok;
	    # canonical definition of a readthrough is
	    # consecutively-ordered genes on the same strand.
	    #
	    # TO DO: also check for intervening genes!  Can still be
	    # readthrough if these are low-expressed, a nice annotation
	    # to have in any case

	    my $distance;
	    # distance from end of geneA to start of geneB.
	    # If geneB comes after geneA on the same strand, should
	    # be a positive number.
	    if (not($pair_strand_ok)) {
	      # genes on different strands, ugh
	      # canonical ordering check is meaningless in this case

	      my @genomic_order = sort {$a->{txStart} <=> $b->{txStart}} ($ra, $rb);
	      dump_die($genomic_order[0], "first", 1);
	      dump_die($genomic_order[1], "second", 1);
	      my $end_first = $genomic_order[0]->{txEnd} || die;
	      my $start_second = $genomic_order[1]->{txStart} || die;

	      if ($start_second > $end_first) {
		# genes are separate under genomic ordering
		$distance = $start_second - $end_first;
		# hack/placeholder as gene model orientations not consistent
	      } else {
		$distance = 0;
		$has_gene_overlap = 1;
		# some genes may overlap, e.g. ADAMDEC1-LOC101929294-01
	      }

	    } elsif ($ra->{strand} eq "+") {
	      # both genes on +
	      my $end_A = $ra->{txEnd};
	      my $start_B = $rb->{txStart};

	      # TO DO: interbase correction, but close enough

#	      printf STDERR "%s %s=%d-%d %s=%d-%d\n",
#		$ra->{chrom},
#		  $geneA_rf, $ra->{txStart}, $ra->{txEnd},
#		    $geneB_rf, $rb->{txStart}, $rb->{txEnd};

	      $distance = $start_B - $end_A;
 	    } elsif ($ra->{strand} eq "-") {
	      # both genes on -
	      my $end_A = $ra->{txStart};
	      my $start_B = $rb->{txEnd};
	      $distance = $end_A - $start_B;
	      # UCSC always orders gene mappings genomically,
	      # so txStart/txEnd are reversed for - strand genes

#	      printf STDERR "%s %s=%d-%d %s=%d-%d dist=%d\n",
#		$ra->{chrom},
#		$geneA_rf, $ra->{txStart}, $ra->{txEnd},
#		  $geneB_rf, $rb->{txStart}, $rb->{txEnd},
#		    $distance;
 	    } else {
	      die;
	    }

	    if ($verbose) {
	      printf STDERR "%s %s: %d\n", $geneA_raw, $geneB_raw, abs($distance);
	    }

	    my $distance_ok;
	    if (abs($distance) <= $max_distance_near) {
	      $has_distance_near = $distance_ok = 1;
	    } elsif (abs($distance) <= $max_distance_far) {
	      $has_distance_far = $distance_ok = 1;
	    }

	    if ($distance_ok) {
	      #
	      # genes are close enough to be potential readthrough
	      #
	      $has_a_b_ordering = 1 if $distance >= 0;
	      # canonical behavior, gene A ordered before gene B on strand.
	      # exceptions may still be of interest however.  Examples:
	      #   AFMID_TNRC6C (+)
	      #   ABR_YWHAE (-)

	      push @rt_rf_pairs, [ $ra, $rb ];
	    }
	  } elsif ($verbose) {
	    printf STDERR "%s %s: chrom mismatch\n", $geneA_raw, $geneB_raw;
	  }


	}  # $rb
      }  #  $ra
    } else {
      printf STDERR "ERROR: can't find refFlat entry for %s\n", $geneA_raw unless $geneA_rf;
      printf STDERR "ERROR: can't find refFlat entry for %s\n", $geneB_raw unless $geneB_rf;
    }
  }

  #
  #  genomic interval spanning all isoform mappings qualifying for RT:
  #
  my $rt_interval = "";
  if (@rt_rf_pairs) {
    my %c2i;
    # save just one interval per chrom
    foreach my $pair (@rt_rf_pairs) {
      my %pos;
      my $chr;
      foreach my $r (@{$pair}) {
	$chr = $r->{chrom} || die;
	# same for all
	foreach my $f (qw(txStart txEnd)) {
	  $pos{$r->{$f}} = 1;
	}
      }
      my @sorted = sort {$a <=> $b} keys %pos;
      my $interval = sprintf '%s:%d-%d', $chr, $sorted[0], $sorted[$#sorted];
      $c2i{$chr} = $interval;
    }

    my @i = map {$c2i{$_}} sort keys %c2i;
    $rt_interval = join ",", @i;
  }
  $self->rt_interval($rt_interval);


  if (($has_distance_near or $has_distance_far) and $contig) {
    # - blat contig vs genome
    # - iterate through each qualifying refflat row pair for A/B
    # - check qualifications for each pair
    # - stop if we find one

    #
    #  get set of chromosomes the genes are mapped to in refFlat:
    #
    my %usable_chrom;
    foreach my $rf_pair (@rt_rf_pairs) {
      my ($rf_a, $rf_b) = @{$rf_pair};
      foreach my $rf (@{$rf_pair}) {
	my $chr = cook_chromosome_name($rf->{chrom} || die);
	$usable_chrom{$chr} = 1;
      }
    }

    my $parser;
    my $tfw = new TemporaryFileWrangler();
    if (0) {
      # debug only
      print STDERR "DEBUG: using test blat result\n";
      $parser = Bio::SearchIO->new(
				   -file => "ARRDC1_EHMT1.blat",
				   -format => "psl"
				  );
    } else {
      my $f_blat_query = $tfw->get_tempfile("-append" => ".fa");
      my $f_blat_out = $tfw->get_tempfile("-append" => ".psl");
      unlink($f_blat_query, $f_blat_out);

      open(BTMP, ">" . $f_blat_query) || die;
      printf BTMP ">contig\n%s\n", $contig;
      close BTMP;

      if ($verbose) {
	my $f_copy = "blat_query.fa";
	printf STDERR "blat query copied to %s\n", $f_copy;
	copy($f_blat_query, $f_copy) || die;
      }

      $bc->outfile($f_blat_out);

      $parser = $bc->query("-fasta" => $f_blat_query);
    }

    my $result = $parser->next_result();
    # one result per query sequence, in this case we just have the one

    my $hit;
    while (my $h = $result->next_hit()) {
      my $hit_name = cook_chromosome_name($h->name());
      my $is_usable = $usable_chrom{$hit_name} || 0;
#      print STDERR "hit $h $hit_name $is_usable\n";
      if ($is_usable) {
	# use first hit to a chrom the gene is mapped to
	$hit = $h;
	last;
      }
    }

    if ($hit) {
      #
      # usable hit
      #

      my $hit_name = $hit->name();
      # chromosome matched
      my $hsp = $hit->next_hsp();
      # for this use case we expect a single HSP to span both genes.

      printf STDERR "   name:%s score:%s strand:%s q:%d-%d subj:%d-%d num_identical:%d frac_identical_query:%s query_span:%d ref_span:%d total_span=%d query_string=%s hit_string=%s\n",
	$hit_name,
	  $hsp->score,
	    $hsp->strand,
	      $hsp->range("query"),
		$hsp->range("hit"),
		  $hsp->num_identical(),
		    $hsp->frac_identical("query"),
		      $hsp->length("query"),
			$hsp->length("hit"),
			  $hsp->length("total"),
			    $hsp->query_string(),
			      $hsp->hit_string() if $verbose;

      #
      #  build a list of block positions the blat hit covers in the target:
      #
      my @block_pos;
      printf STDERR "gap blocks:\n" if $verbose;
      foreach my $k (qw(query hit)) {
	printf STDERR "  %s\n", $k if $verbose;
	foreach my $gb (@{$hsp->gap_blocks($k)}) {
	  printf STDERR "    %s\n", join ",", @{$gb} if $verbose;

	  if ($k eq "hit") {
	    my ($start, $len) = @{$gb};
	    my $pos = $start + int($len / 2);
	    # pick a point solidly in the center of the hit, in case
	    # there's a little microhomology, etc. at the edge which
	    # could potentially interfere with matches to exons
	    push @block_pos, $pos;
	  }

	}
      }

      if (@block_pos) {
	# span found

	#
	#  annotate gene and exon # for each aligned block:
	#
	foreach my $rf_pair (@rt_rf_pairs) {
	  # foreach refFlat record pair for genes A and B...
	  my ($rf_a, $rf_b) = @{$rf_pair};

	  my $strand = $rf_a->{strand} || die;

	  my @block2gene;
	  # block matches, ordered GENOMICALLY rather than by strand

	  for (my $i = 0; $i < @block_pos; $i++) {
	    my ($type_a, $fno_a, $strand_a) = $rf->get_annotation_for_position(
									       "-row" => $rf_a,
									       "-base" => $block_pos[$i],
									       "-intergenic-ok" => 1,
									       "-extended" => 1
									      );
	    my ($type_b, $fno_b, $strand_b) = $rf->get_annotation_for_position(
									       "-row" => $rf_b,
									       "-base" => $block_pos[$i],
									       "-intergenic-ok" => 1,
									       "-extended" => 1
									      );

	    if ($type_a eq "coding") {
	      $block2gene[$i] = [ "A", $fno_a ];
	    } elsif ($type_b eq "coding") {
	      $block2gene[$i] = [ "B", $fno_b ];
	    }
	  }

	  #
	  #  find the boundary where the block gene annotation changes:
	  #
	  for (my $i = 0; $i < @block2gene - 1; $i++) {
	    if ($block2gene[$i] and $block2gene[$i + 1] and
		# we know what gene belongs to each block
		$block2gene[$i]->[0] ne $block2gene[$i + 1]->[0]) {

	      my $first_is_2nd_to_last_exon;
	      my $second_is_exon_2;
	      my $ecount_a = $rf_a->{exonCount} || die;

	      if ($strand eq "+") {
		# genomic block order matches gene order
		$first_is_2nd_to_last_exon = $block2gene[$i]->[1] == $ecount_a - 1;
		$second_is_exon_2 = $block2gene[$i + 1]->[1] == 2;
	      } elsif ($strand eq "-") {
		$first_is_2nd_to_last_exon = $block2gene[$i + 1]->[1] == $ecount_a - 1;
		$second_is_exon_2 = $block2gene[$i]->[1] == 2;
	      }

	      if ($first_is_2nd_to_last_exon and $second_is_exon_2) {
		$has_second_to_last_to_second = 1;
	      }
	    }
	  }
	}			# $rf_pair

      }
    }
  }

  my $is_rt = "";
  my @rt_info;
  if ($has_distance_near or $has_distance_far) {
    if ($has_distance_near and $has_consistent_strand and $has_a_b_ordering) {
      # meets canonical definition
      $is_rt = "yes";
    } else {
      # possible, but some details don't appear canonical
      $is_rt = "maybe";
    }

    push @rt_info, "second_to_last_exon_to_second_exon" if $has_second_to_last_to_second;
    push @rt_info, "distant" unless $has_distance_near;
    push @rt_info, "inconsistent_strand" unless $has_consistent_strand;
    push @rt_info, "gene_mappings_overlap" if $has_gene_overlap;
    push @rt_info, "inconsistent_gene_ordering" unless $has_a_b_ordering;
  }
  $self->rt_info(\@rt_info);

  return $is_rt;
}

1;

# LRF support:
# _______________      ______________________________      _______________
#                \____/                              \____/
