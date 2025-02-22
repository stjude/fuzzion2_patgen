package GenomeUtils;
# MNE 8/2013

use strict;
use Carp qw(confess);
use File::Copy;
use Cwd qw(abs_path);

use FileUtils qw(universal_open find_binary);
use SampleTumorNormal;
use FAI;
use MiscUtils qw(dump_die);

@GenomeUtils::ISA = qw(Exporter);
@GenomeUtils::EXPORT_OK = qw(
                           cook_chromosome_name
reverse_complement
complement

check_reference_name_compatibility
parse_interval
bam_viewer_link
merge_sorted_bams

convert_bambino_insertion_position
                          );

%GenomeUtils::NC2CHR = (
	      #
	      # human reference assembly:
	      #
	      "NC_000001" => 1,
	      "NC_000002" => 2,
	      "NC_000003" => 3,
	      "NC_000004" => 4,
	      "NC_000005" => 5,
	      "NC_000006" => 6,
	      "NC_000007" => 7,
	      "NC_000008" => 8,
	      "NC_000009" => 9,
	      "NC_000010" => 10,
	      "NC_000011" => 11,
	      "NC_000012" => 12,
	      "NC_000013" => 13,
	      "NC_000014" => 14,
	      "NC_000015" => 15,
	      "NC_000016" => 16,
	      "NC_000017" => 17,
	      "NC_000018" => 18,
	      "NC_000019" => 19,
	      "NC_000020" => 20,
	      "NC_000021" => 21,
	      "NC_000022" => 22,
	      "NC_000023" => "X",
	      "NC_000024" => "Y",
              # these have historically been stable (GRCh37)

              "NC_060925" => 1,
               # hs1

	      "NC_001807" => "MT",
	      "NC_012920" => "MT",   # 03-APR-2023
              # ...however mitochondrion accession may vary
    );

my $VERBOSE_NAME_MUNGE = 0;

sub cook_chromosome_name {
  my ($chr, %options) = @_;

  my $ucsc = $options{"-ucsc"};
  my $genome = $options{"-genome"};

  my ($munge_23, $munge_24, $munge_25);
  if ($genome) {
#    if ($genome =~ /^hg/ or $genome eq "GRCh37-lite") {
    if ($genome =~ /^hg/ or $genome =~ /GRCh\d+/) {
      # human
      $munge_23 = $munge_24 = $munge_25 = 1;
    } elsif ($genome =~ /^Zv/i) {
      # zebrafish
      $munge_23 = $munge_24 = $munge_25 = 0;
    } elsif (
	     $genome =~ /^mm/i or
	     $genome =~ /^dm/i or
	     $genome eq "MGSCv37" or
	     $genome eq "BDGPr5"
	    ) {
      # TEMPORARY HACK:
      # needs to be fixed (fewer chroms for mouse!)
      $munge_23 = $munge_24 = $munge_25 = 1;
    } else {
      confess "unhandled genome $genome";
    }
  } else {
    # assume human
    $munge_23 = $munge_24 = $munge_25 = 1;
  }


  my $chr_raw = $chr;
  $chr =~ s/^chr//i;
  if ($chr =~ /^\d+$/ and $chr >= 1 and $chr <= 22) {
    # ok
    # constants: HACK
  } elsif ($chr eq "23") {
    # conversion of string chromosome names to integer?
    # observed in COSMIC
    if ($munge_23) {
      $chr = "X";
      printf STDERR "WARNING: mapping %s => %s\n", $chr_raw, $chr if $VERBOSE_NAME_MUNGE;
    }
  } elsif ($chr eq "24") {
    # conversion of string chromosome names to integer?
    # observed in COSMIC
    if ($munge_24) {
      $chr = "Y";
      printf STDERR "WARNING: mapping %s => %s\n", $chr_raw, $chr if $VERBOSE_NAME_MUNGE;
    }
  } elsif ($chr eq "25") {
    # conversion of string chromosome names to integer?
    # observed in COSMIC
    if ($munge_25) {
      $chr = "MT";
      printf STDERR "WARNING: mapping %s => %s\n", $chr_raw, $chr;
    }
  } elsif ($chr eq "X" or $chr eq "Y") {
    # ok
  } elsif ($chr eq "M" or $chr eq "MT") {
    $chr = "MT";
  } elsif ($chr =~ /^NC_/) {
    my ($acc, $version) = split /\./, $chr;
    $chr = $GenomeUtils::NC2CHR{$acc} || die "can't map $acc to chrom";
    # consider this a fatal error: while primary chromosome accessions
    # are stable, mitochondrion may require updates (see NC2CHR above)
  } elsif ($options{"-return-unknown"}) {
    # if not in known format, silently return original
    $chr = $chr_raw;
  } elsif ($options{"-undef-unknown"}) {
    $chr = undef;
  } else {
    confess "unknown chrom $chr";
  }

  if ($ucsc) {
    $chr = "M" if $chr eq "MT";
    $chr = "chr" . $chr;
  }

  return $chr;
}

sub reverse_complement {
  my $source = shift;
  if (ref $source eq "ARRAY") {
    my @r = reverse @{$source};
    foreach (@r) {
#      tr/ACGTacgt/TGCAtgca/;
#      tr/acgtACGTrymkRYMK/tgcaTGCAyrkmYRKM/;
      # S and W remain the same

      tr/acgtACGTrymkRYMKBVDHbvdh/tgcaTGCAyrkmYRKMVBHDvbhd/;
    }
    return @r;
  } else {
#    (my $result = reverse $source) =~ tr/ACGTacgt/TGCAtgca/;
#    (my $result = reverse $source) =~ tr/acgtACGTrymkRYMK/tgcaTGCAyrkmYRKM/;
    (my $result = reverse $source) =~ tr/acgtACGTrymkRYMKBVDHbvdh/tgcaTGCAyrkmYRKMVBHDvbhd/;
    return $result;
  }
}

sub complement {
  my ($thing) = @_;
  if (ref $thing) {
    die "not implemented";
  } else {
    $thing =~ tr/acgtACGTrymkRYMKBVDHbvdh/tgcaTGCAyrkmYRKMVBHDvbhd/;
  }
  return $thing;
}

sub check_reference_name_compatibility {
  my (%options) = @_;
  # TO DO: extend for other types
  my $fasta = $options{"-fasta"};
  my $vcf = $options{"-vcf"};

  my $set1;
  my $set2;

  if (0) {
    print STDERR "DEBUG!!\n";
    $set1 = [ 1 .. 22, qw(X Y MT)];
    $set2 = [ 1 .. 22, qw(X Y M)];
  } elsif ($fasta and $vcf) {
    my $fai = $fasta . ".fai";
    die "where is $fai" unless -s $fai;
    $set1 = [];
    open(RNCTMP, $fai) || die;
    while (<RNCTMP>) {
      chomp;
      my @f = split /\t/, $_;
      push @{$set1}, $f[0];
    }

    my $fh_tmp = universal_open($vcf) || die;
    my %v;
    printf STDERR "scanning %s...", $vcf;
    while (<$fh_tmp>) {
      next if /^\#/;
      my $ref = (split /\t/, $_)[0];
      printf STDERR "%s...", $ref unless $v{$ref};
      $v{$ref} = 1;
#      $v{$f[0]} = 1;
    }
    print STDERR "\n";
    $set2 = [sort keys %v];
  }

  my %trouble;
  foreach my $set ($set1, $set2) {
    foreach my $ref (@{$set}) {
      my $cooked = cook_chromosome_name($ref, "-return-unknown" => 1);
      $trouble{$cooked}{$ref} = 1;
    }
  }

  my $ok = 1;
  foreach my $cooked (sort keys %trouble) {
    my $raw = $trouble{$cooked};
    if (scalar keys %{$raw} > 1) {
      printf STDERR "ERROR for %s: observed %s\n", $cooked, join ",", sort keys %{$raw};
      $ok = 0;
    }
  }
  return $ok;
}

sub parse_interval {
  my ($interval) = @_;
  my @f = split /:/, $interval;
  my $error_msg = "interval $interval must be in format X:123[-456]";
  die $error_msg unless @f == 2;
  my ($chr, $region) = @f;
  @f = split /\-/, $region;
  my ($start, $end);
  if (@f == 2) {
    # range
    ($start, $end) = @f;
  } elsif (@f == 1) {
    # single position
    $start = $end = $f[0];
  } else {
    die $error_msg;
  }
  foreach ($start, $end) {
    confess "position $_ in $interval not numeric" unless /^\d+$/;
  }

  if (wantarray) {
    return ($chr, $start, $end);
  } else {
    die "parse_interval: interval scalar/hash result";
  }
}

sub bam_viewer_link {
  my %options = @_;
  my $bams = $options{"-bams"} || die "-bams";
  my $is_clinical = $options{"-clinical"};
  my $genome = $options{"-genome"} || "hg19";
  if ($genome eq "GRCh37-lite" or $genome eq "hg19") {
    $genome = "hg19";
  } elsif ($genome eq "GRCh38") {
    $genome = "hg20";
    # no such thing, but this is the name used by the server
  }
      

  my %link;
  $link{ref} = $genome;
  $link{region} = $options{"-chr"} || die "-chr";
  $link{center} = $options{"-center"} || die "-center";
  $link{fullPath} = "true";

  my $stn = new SampleTumorNormal();
  my %bams;
  my @bams;
  if (ref $bams eq "ARRAY") {
    @bams = @{$bams};
  } else {
    # single file
    @bams = $bams;
  }
  die "no bams" unless @bams;

  foreach my $bam (@bams) {
    $bam = abs_path($bam);
    # some sample renaming hacks involve temporary symlinks in home
    # directories; resolve to final path so these will work with viewer

    $is_clinical = 1 if $bam =~ /cgs01/ or $bam =~ /clingen/;
    # hack
    if ($stn->parse($bam)) {
      if ($stn->is_tumor) {
	$bams{tumorname} = $bam;
      } elsif ($stn->is_normal) {
	$bams{normalname} = $bam;
      } else {
	die "can't determine tumor/normal for $bam";
      }
    } else {
      printf STDERR "WARNING: can't determine tumor/normal status for %s, assuming T\n", $bam;
      $bams{tumorname} = $bam;
    }
  }

  my $link = sprintf 'http://bamviewer-%s:8080/BAMViewer/aceview/splash?', $is_clinical ? "cp" : "rt";
  die unless %bams;
  foreach (sort keys %bams) {
    $link .= sprintf '%s=%s&', $_, $bams{$_};
  }

  foreach my $k (sort keys %link) {
    $link .= sprintf '%s=%s&', $k, $link{$k};
  }
  $link =~ s/&$//;

  return $link;
  
}

sub merge_sorted_bams {
  my (%options) = @_;
  my $f_out = $options{"-out"} || die "-out";
  my $bams = $options{"-bams"} || die "-bams";
  die unless @{$bams};

  if (@{$bams} == 1) {
    # samtools merge requires at least 2 files
    copy($bams->[0], $f_out) || die "can't copy to $f_out";
  } else {
    my $f_out_tmp = $f_out . ".tmp";
    my $cmd = sprintf 'samtools merge %s %s', $f_out_tmp, join " ", @{$bams};
    printf STDERR "running: %s...\n", $cmd;
    system $cmd;
    die "$cmd exited with $?" if $?;
    die unless -s $f_out_tmp;
    rename($f_out_tmp, $f_out) || die "can't rename $f_out_tmp to $f_out";
  }

  my $cmd = "samtools index $f_out";
  printf STDERR "running: %s...\n", $cmd;
  system $cmd;
  die "$cmd exited with $?" if $?;
}

sub convert_bambino_insertion_position {
  # Bambino insertions use the base number AFTER the inserted sequence.
  # Given variant alleles, if the event is an insertion, adjust
  # to the base number before the insertion.  Helpful for some
  # comparisons and for batch tabix variant databases, which expect
  # 0-based.
  my ($pos, $ref, $alt) = @_;
  foreach ($ref, $alt) {
    $_ = "" if $_ eq "-";
    # SJ post
  }
  $pos-- if length($ref) == 0 and length($alt);
  return $pos;
}

1;
