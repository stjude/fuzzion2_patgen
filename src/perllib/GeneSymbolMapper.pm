package GeneSymbolMapper;
# attempt to translate gene symbols in one set to another
# MNE 4/2014

# tests:
#
# SV:
#
#
# CNV analysis using SV genes in GENE_EXON_REGION:
# MADH4: should be SMAD4 (later symbol)
# WTX:  updated symbol is AMER1 but GENE_EXON_REGION uses FAM123B (older)
#

use strict;
use Carp qw(confess);

use Configurable;
use Exporter;

use MiscUtils qw(dump_die);
use HGNCParser;
use GeneSymbolStandardizer;
use TdtConfig;

@GeneSymbolMapper::ISA = qw(Configurable Exporter);
@GeneSymbolMapper::EXPORT_OK = qw(new_gsm_lite);

use MethodMaker qw(
	genes
	hgnc
	gss

	hgnc_file
	eg_file
	refseq2gene

problem
approved_symbol

hgnc_synonym_prune
hgnc_previous_prune
hgnc_synonym_enable

enable_entrez_gene
hgnc_lite

clone
genes_ordered
user_hash
unique_approved
		  );

sub new {
  my ($type, %options) = @_;
  my $self = {};
  bless $self, $type;
  $self->enable_entrez_gene(1);
  # as of 2018 may no longer be necessary
  $self->hgnc_synonym_enable(1);
  $self->hgnc_synonym_prune(1);
  # enable HGNC synonyms by default, but prune some synonyms
  # to reduce risk of dubious matches.
  #
  # example:
  # - existing germline reviewable list contains FANCA
  # - new reviewable gene FAH
  # - both FAH and FANCA are approved symbols
  # - however, FAH is also be a synonym for FANCA
  # - these are two completely different genes, a synonym match
  #   is not strong enough evidence to reject gene as a duplicate
  $self->hgnc_previous_prune(1);
  # however, some old symbols map to more than one new symbol!
  $self->reset_instance();
  $self->configure(%options);
  if (my $parent = $self->clone()) {
    %{$self} = %{$parent};
    $self->reset_instance();
  } else {
    $self->setup();
  }
  return $self;
}

sub reset_instance {
  # reset user mapping database
  my ($self) = @_;
  $self->genes({});
  $self->refseq2gene({});
  $self->genes_ordered([]);
  $self->user_hash({});
}

sub setup {
  my ($self) = @_;

  my $hgnc = new HGNCParser(
			    ("-file" => $self->hgnc_file || die "-hgnc_file"),
			    "-lite" => $self->hgnc_lite(),
			   );
  $hgnc->prune_synonyms() if $self->hgnc_synonym_prune();
  $hgnc->prune_previous() if $self->hgnc_previous_prune();
  $self->hgnc($hgnc);
  # HGNC/HUGO

  if ($self->enable_entrez_gene) {
    # Entrez Gene
    $self->gss(new GeneSymbolStandardizer("-filename" => $self->eg_file || die "-eg_file"));
  }
}

sub add_gene {
  my ($self, %options) = @_;
  my $gene = $options{"-gene"};
  confess "-gene" unless defined $gene;
  my $refseq = $options{"-refseq"};
  # optional
  my %record;

  if ($gene =~ /^(C\d+)ORF(\d+)$/) {
    # capitalized "ORF" in gene name, e.g. from COSMIC
    my $fixed = $1 . "orf" . $2;
    printf STDERR "WARNING: reformatting %s as %s\n", $gene, $fixed;
    $gene = $fixed;
  }


  $record{gene} = $gene;
  $record{refseq} = $refseq;
  push @{$self->genes_ordered()}, $gene unless $self->genes()->{$gene};
  # only add to ordered list once

  if ($refseq) {
    my $r2g = $self->refseq2gene();
    if (my $existing = $r2g->{$refseq}) {
      printf STDERR "WARNING: duplicate refGene entry for %s: old:%s new:%s\n", $refseq, $existing, $gene unless $existing eq $gene;
      # can happen for hack NG_ entries, see COMPBIO-2772
    } else {
      $r2g->{$refseq} = $gene;
    }
  }

  push @{$self->genes->{$gene}}, \%record;
}


sub reset {
  my ($self) = @_;
  $self->problem("");
  $self->approved_symbol("");
  # reset before sift_hgnc_hits() called for the first time
}

sub resolve {
  my ($self, %options) = @_;
  my $sym = $options{"-symbol"};
  confess "-symbol" unless defined $sym;
  my $result;
  $self->reset();

  if (0 and $self->contains($sym)) {
    # target is the same in both lists.
    # normally we'd check before calling this routine and this would
    # be a fatal error as a lookup is not needed, and we don't
    # want to require database matching/filtering if so.
    #
    # However for some use cases (e.g. universal somatic gene list
    # construction) it's desirable to run these symbols through
    # to determine the HUGO symbol for each, i.e. to detect
    # synonymous symbols in combined list.
    die "same";
  }

  my $hgnc = $self->hgnc();
  my $problem;
  my $method;

  if (0) {
    my $test = $hgnc->find("-symbol" => $sym,
			   "-approved" => 1);
    dump_die($test->[0]);
  }

  unless ($result or $problem) {
    ($result, $problem) = $self->sift_hgnc_hits(%options,
						"-label" => "approved $sym",
						"-hits" =>
						$hgnc->find("-symbol" => $sym,
							    "-approved" => 1)
					       );
    $method = "approved" if $result;
  }
  # ACKR3: a SV gene not found in GENE_EXON_REGION,
  # but NM_020311 and older symbol CXCR7 are

  unless ($result or $problem) {
    ($result, $problem) = $self->sift_hgnc_hits(%options,
						"-label" => "previous $sym",
						"-hits" =>
						$hgnc->find("-symbol" => $sym,
							    "-previous" => 1)
					       );
    $method = "previous" if $result;
    # CEP1 to GENE_EXON_REGION via NM_007018:
    # CEP1 -> CNTRL -> NM_007018 -> CEP110
    # this one is tricky because CEP1 is also a synonym for CDC42EP1.
    # avoid by using previous symbol lookups before synonym lookups.
  }

  if ($self->hgnc_synonym_enable and not($result or $problem)) {
    ($result, $problem) = $self->sift_hgnc_hits(%options,
						"-label" => "synonym $sym",
						"-hits" =>
						$hgnc->find("-symbol" => $sym,
							    "-synonym" => 1)
					       );
    $method = "synonym" if $result;
  }

  unless ($result or $problem) {
    # e.g. attempting to map NPSR1-AS1, a new symbol, to AAA1 in 
    # GENE_EXON_REGION.  Lookup works in reverse, i.e. AAA1 is a synonym
    # for NPSR-AS1l, but NPSR-AS11 doesn't work in any index to AAA1.
    # So, if all else fails, and we can find a HGNC record for the lookup
    # symbol, search for previous and synonym symbols.

    if (my $r = $hgnc->find_hgnc_record("-symbol" => $sym)) {
      my @prev = split /,\s*/, $r->{"Previous Symbols"};
      my @syn = split /,\s*/, $r->{"Synonyms"};

      foreach my $sym (@prev) {
	next if $hgnc->is_blacklisted_previous($sym);
	if ($self->contains($sym)) {
	  $result = $sym;
	  $method = "rescue_previous";
	  last;
	}
      }
      foreach my $sym (@syn) {
	next if $hgnc->is_blacklisted_synonym($sym);
	if ($self->contains($sym)) {
	  $result = $sym;
	  $method = "rescue_synonym";
	  last;
	}
      }
    }
  }

  if ($self->enable_entrez_gene()) {
    unless ($result or $problem) {
      # try Entrez Gene if all else fails, occasionally has
      # an entry (e.g. ALO17 -> RNF213)
      my $syns = $self->gss->find($sym);
      if ($syns and ref $syns) {
	if (@{$syns} == 1 and $self->contains($syns->[0])) {
	  ($result) = @{$syns};
	  $method = "entrez_gene";
	} else {
	  $problem = 1;
	}
      }
    }
  }

#  printf STDERR "no luck for %s\n", $sym unless $result;


#     if ($hugo_hits) {
#       if (@{$hugo_hits} == 1) {
# 	# only 1 result
# 	my $approved = $hugo_hits->[0]->{"Approved Symbol"} || die;
# 	if ($self->contains($approved)){
# 	  # given symbol translates to an approved symbol which
# 	  # is found in the match list
# 	  $result = $approved;
# 	} else {
# 	  die "can't find approved $approved";
# 	  # try older symbols??
# 	}
#       } else {
# 	# ambiguous matches, e.g.
# 	# HSPCA, which hits HSP90AA1 and HSP90AA2.
# 	my @try = map {$_->{"Approved Symbol"}} @{$hugo_hits};
# #	die @try;
# 	my @matches;
# 	foreach my $g (@try) {
# 	  push @matches, $g if $self->contains($g);
# 	}

# 	if (@matches == 1) {
# 	  ($result) = @matches;
# 	} else {
# 	  die "not matches or ambiguous";
# 	}
#       }
#     }

  if ($result and $ENV{GSM_VERBOSE}) {
    printf STDERR "GSM lookup %s: hit %s via %s\n", $sym, $result, $method;
  }

  return $result;
}

sub contains {
  my ($self, $sym) = @_;
  return $self->genes->{$sym} ? 1 : undef;
}

sub sift_hgnc_hits {
  my ($self, %options) = @_;
  die unless exists $options{"-hits"};
  my $hugo_hits = $options{"-hits"};
  my $hgnc = $self->hgnc();
  my $result;
  my %found;
  my $label = $options{"-label"} || "";
  if ($hugo_hits) {

    my %approved;

    foreach my $hh (@{$hugo_hits}) {
      my $approved = $hh->{"Approved Symbol"} || die;
      $approved{$approved} = 1;

      my $match;
      if ($self->contains($approved)) {
	# latest approved symbol is in the target set
	$found{$approved} = 1;
	$match = 1;
      }

      if (not($match) and $hh->{refseqs_hash}) {
	# the latest approved symbol isn't in the target set,
	# however it does contain NM records, try that lookup.
	foreach my $rs (keys %{$hh->{refseqs_hash}}) {
	  my $hit = $self->refseq2gene->{$rs};
	  if ($hit) {
	    $found{$hit} = 1;
	    $match = 1;
	  }
	}
      }

      if (not($match) and my $ps = $hh->{"Previous Symbols"}) {
	# see if any of the previous symbols for this entry 
	# are found in the target set, e.g.
	#   1. source symbol (SV config): CMKOR1 ->
	#   2. approved symbol (FB): ACKR3 ->
	#   3. target symbol (GER): CXCR7
	my @prev = split /,\s*/, $ps;
	foreach my $sym (@prev) {
	  if (my $reason = $hgnc->is_blacklisted_previous($sym)) {
	    # e.g. TCF3 / TXF7L1: if searching for TXF7L1 while the TCF3
	    # entry will have been pruned from the previous symbol index,
	    # it will pop up again during this check, but should be ignored.
#	    printf STDERR "ignoring blacklisted previous symbol %s (%s)\n", $sym, $reason;
	  } elsif ($self->contains($sym)) {
	    $found{$sym} = 1;
	    $match = 1;
	  }
	}
      }
    }

    if (%approved and !$self->approved_symbol()) {
      my $list = join ",", sort keys %approved;
      $self->approved_symbol($list);
      # only record first/best match
      # there may be more than one, e.g. lookup of SIL:
      # - match found via synonyms only
      # - hits multiple genes: STIL, PMEL
    }

  }

  my $problem;

  if (my @found = keys %found) {
    if (@found == 1) {
      ($result) = keys %found;
#      printf STDERR "found %s => %s\n", $label, $result if $label and 1;
    } else {
      $problem = "problematic_ambiguity";
      printf STDERR "ERROR: %s for %s: %s\n", $problem, $options{"-symbol"}, join ",", @found;
    }
  }

  return ($result, $problem);
}

sub populate_refflat {
  my ($self, %options) = @_;
  my $rf = $options{"-refflat"} || die "-refflat";

  open(RFTMP, $rf) || die;
  my $genes = 0;
  while (<RFTMP>) {
    chomp;
    my @f = split /\t/, $_;
    my $gene = clean_sharp_gene_symbol($f[0]);
    if ($gene) {
      $genes++;
      my $nm = $f[1];
      $self->add_gene(
		     "-gene" => $gene,
		     "-refseq" => $nm
		    );
    } else {
#      printf STDERR "no gene symbol in refFlat entry %s\n", $_;
      # e.g. _locPar NR_001526-locPar
    }
  }
  close RFTMP;
  die "no genes" unless $genes;
}


sub clean_sharp_gene_symbol {
  # in Michael Rusch's "sharp" versions of refFlat,
  # ambiguous loci assigned "_loc" suffix, e.g. DUX4_locF.
  # These will also be reported by FusionBuilder with this suffix,
  # so need to remove it to get the raw gene symbol.
  my ($gene) = @_;
  if ($gene =~ /_loc/) {
    my $before = $gene;
    $gene =~ s/_loc.*$//;
#    printf STDERR "stripping %s => %s\n", $before, $gene;
  }
  return $gene;
}

sub populate_gene_exon_region {
  #
  # symbols reported by GENE_EXON_REGION
  #
  my ($self, %options) = @_;
  my $ger_dir = $options{"-dir"} || die;
  my @files = glob("$ger_dir/*txt");
  die unless @files;
  foreach my $f (@files) {
    open(GERTMP, $f) || die;
    while (<GERTMP>) {
      my @f = split /\t/, $_;
      my ($gene, $nm) = @f[0,2];
      $self->add_gene(
		     "-gene" => $gene,
		     "-refseq" => $nm
		    );
    }
  }
  close GERTMP;
}

sub populate_chr2gene {
  my ($self, %options) = @_;
  my $dir = $options{"-dir"} || die;
  my @files = glob($dir . "/*gene.txt");
  die "no files in $dir" unless @files;
  my %saw;
  foreach my $f (@files) {
    open(CGTMP, $f) || die;
    while(<CGTMP>) {
      chomp;
      my @f = split /\t/, $_;
      die unless @f == 4;
      my @list = split /\|/, $f[0];
      my ($gene, $nm) = @list;
      unless ($saw{$gene}) {
	$self->add_gene(
			"-gene" => $gene,
			"-refseq" => $nm
		       );
	$saw{$gene} = 1;
      }
    }
  }
}

sub find {
  # highest-level search
  my ($self, $sym) = @_;
  my $result;
  if ($self->contains($sym)) {
    $result = $sym;
  } else {
    $result = $self->resolve("-symbol" => $sym);
  }
  return $result;
}

sub new_gsm_lite {
  # instance factory from clone
  my (%options) = @_;

  unless ($GeneSymbolMapper::GSM_MAIN) {
    my $lite_mode = exists $options{"-lite"} ? $options{"-lite"} : 1;
    my $f_hgnc;
    if ($f_hgnc = $options{"-hgnc"} || $ENV{GSM_HGNC}) {
      printf STDERR "using %s\n", $f_hgnc;
    } else {
      my $config_species = TdtConfig::readConfig('species', "Homo_sapiens") || die "can't find hs species config";
      $f_hgnc = $config_species->{HGNC} || die "no HGNC config";
    }

    $GeneSymbolMapper::GSM_MAIN  = new GeneSymbolMapper(
      "-hgnc_file" => $f_hgnc,
      "-hgnc_lite" => $lite_mode,
      "-enable_entrez_gene" => 0,
	);
  }

  return new GeneSymbolMapper("-clone" => $GeneSymbolMapper::GSM_MAIN);
}

sub get_hgnc_approved {
  my ($self) = @_;
  my @genes = grep {!/~withdrawn$/i} keys %{$self->hgnc->index_approved()};
  return \@genes;
}

sub get_gene_count {
  my ($self) = @_;
  return scalar keys %{$self->genes()};
}

sub get_gene_list {
  my ($self) = @_;
#  return [ keys %{$self->genes()} ];
  return $self->genes_ordered();
}

sub validate_genes {
  my ($self, %options) = @_;
  my $label = $options{"-label"};
  my $genes = $self->genes_ordered();
  my $hgnc = $self->hgnc();
  my %counts;
  my $verbose = 1;

  my %approved;
  $self->unique_approved(\%approved);

  foreach my $gene (@{$genes}) {
    my $key;
    my $approved;
    if ($hgnc->is_approved($gene)) {
      $key = "approved";
      $approved = $gene;
    } elsif ($hgnc->is_previous($gene)) {
      $key = "previous";
      my $hits = $hgnc->find("-symbol" => $gene,
			    "-previous" => 1);
      $approved = $hits->[0]->{"Approved Symbol"} || die;
    } elsif ($hgnc->is_synonym($gene)) {
      $key = "synonym";
      my $hits = $hgnc->find("-symbol" => $gene,
			    "-synonym" => 1);
      $approved = $hits->[0]->{"Approved Symbol"} || die;
    } else {
      $key = "unknown";
    }
    $counts{$key}{$gene} = 1;
    $approved{$approved} = 1 if $approved;
  }

  if ($verbose) {
    printf STDERR "%s ", $label if $label;
    printf STDERR "gene symbol summary:\n";
    foreach my $key (qw(approved previous synonym unknown)) {
      my $count = scalar keys %{$counts{$key}};
      if ($count) {
	printf STDERR " %10s: %d", $key, $count;
	printf STDERR " (%s)", join ",", sort keys %{$counts{$key}} if $count <= 10;
	print STDERR "\n";
      }
    }
  }

  return $counts{unknown} ? 0 : 1;
}

sub add_user_hash {
  # track a user hash containing a gene symbol field
  my ($self, %options) = @_;
  my $row = $options{"-row"} || die "-row";
  my $f_gene = $options{"-f-gene"} || die "-f-gene";
  my $unique = $options{"-unique"};
  die "-unique" unless defined $unique;
  # require unique key or not

  my $user_hash = $self->user_hash() || die;
  my $gene = $row->{$f_gene} || dump_die($row, "no $f_gene");
#  confess "$gene not unique" if $unique and $self->find($gene);
  confess "$gene not unique" if $unique and $self->contains($gene);
  # require RAW symbols to be unique even if COOKED symbols aren't

  $self->add_gene("-gene" => $gene);
  $user_hash->{$gene} = $row;
}

sub find_user_hash {
  my ($self, $gene_raw) = @_;
  my $result;
  if (my $gene_local = $self->find($gene_raw)) {
#    confess "hey now $gene_local $gene_raw" if $gene_local ne $gene_raw;
    $result = $self->user_hash()->{$gene_local};
  }
  return $result;
}

1;

# LRF support:
# _______________      ______________________________      _______________
#                \____/                              \____/
