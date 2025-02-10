package GeneListMatcher;
# find matches between lists of genes (e.g. fusion events), compensating
# for gene symbol changes
# MNE 11/2021
#
# TO DO:
# - must the genes match in order?  or is reciprocal OK too?

use strict;
use Exporter;

use Configurable;
use MiscUtils qw(get_hash_option dump_die);
use GeneSymbolMapper qw(new_gsm_lite);

@GeneListMatcher::ISA = qw(Configurable Exporter);
@GeneListMatcher::EXPORT_OK = qw();

use MethodMaker qw(
		    gsm
		    raw2set
		    allow_reciprocal
		    match_stats
		  );

sub new {
  my ($type, %options) = @_;
  my $self = {};
  bless $self, $type;
  $self->raw2set({});
  $self->gsm(new_gsm_lite());
  $self->match_stats({});
  $self->configure(%options);
  return $self;
}

sub add_set {
  # add a new event (e.g. fusion)
  my ($self, %options) = @_;
  my $genes = $options{"-genes"} || die "-genes";
  # TO DO: object rather than gene list?
  my %record;
  $record{genes} = [ @{$genes} ];
  $record{data} = $options{"-data"};
  # optional associated data / row

  my $raw2set = $self->raw2set || die;
  my $gsm = $self->gsm || die;
  foreach my $gene (@{$genes}) {
    $gsm->add_gene("-gene" => $gene);
    # disambiguation lookup
    push @{$raw2set->{$gene}}, \%record;
  }
}

sub find {
  my ($self, $genes_query) = @_;
  # find set of all possible matches

  my $gene_count = scalar @{$genes_query};

  my $gsm = $self->gsm || die;
  my $raw2set = $self->raw2set || die;
  my $allow_reciprocal = $self->allow_reciprocal();
  die "allow_reciprocal not set [0|1]" unless defined $allow_reciprocal;

  my %saw;
  my @try;
  # set of lists containing a match to one or more genes in query
  foreach my $gene (@{$genes_query}) {
    if (my $gene_local = $gsm->find($gene)) {
      foreach my $set (@{$raw2set->{$gene_local}}) {
	my $set_gene_count = scalar @{$set->{genes}};

	if ($set_gene_count == $gene_count and
	    # only consider sets with the same number of genes
	    not($saw{$set})
	    # only consider each set once
	   ) {
	  push @try, $set;
	  $saw{$set} = 1;
	}
      }
    }
  }


  my @q_search = $genes_query;
  push @q_search, [ reverse @{$genes_query} ] if $allow_reciprocal;

  #
  # check each candidate set: does EVERY gene in the query match?
  #
  my @hits;

  my $match_stats = $self->match_stats() || die;

  foreach my $set_obj (@try) {
    #
    # check each putative target set
    #
    my $is_reciprocal = 0;
    my $set = $set_obj->{genes} || die;

    foreach my $qs (@q_search) {
      #
      #  vs. query set (also reciprocal, if desired)
      #
      my $result = 1;
      my $perfect = 1;

      for (my $i = 0; $i < $gene_count; $i++) {
	my $gene_q = $qs->[$i];
	my $gene_s = $set->[$i];
	if ($gene_q eq $gene_s) {
	  # identical, continue
	  # TO DO: capitalization?
	} else {
	  my $gsm_tmp = new_gsm_lite();
	  $gsm_tmp->add_gene("-gene" => $gene_s);
	  if (my $alt = $gsm_tmp->find($gene_q)) {
	    # found match via alternate symbol, continue
#	    printf STDERR "$gene_q matches $gene_s via $alt\n";
	    $perfect = 0;
	  } else {
	    # match failure, stop
	    $result = 0;
	    last;
	  }
	}
      }

      if ($result) {
	# match
	push @hits, $set_obj;

	my $type_reciprocal = $is_reciprocal ? "reciprocal" : "primary";
	my $type_gene = $perfect ? "exact" : "alias";

	$set_obj->{match_type_reciprocal} = $type_reciprocal;
	$set_obj->{match_type_gene} = $type_gene;
	$set_obj->{match_genes} = join ",", @{$set};

	my @stat_key;
	push @stat_key, $type_reciprocal;
	push @stat_key, $type_gene;
	my $key = join "_", @stat_key;
	$match_stats->{$key}++;

#	printf STDERR "q:%s target:%s type:%s\n", join(",", @{$genes_query}), join(",", @{$set}), $key;

	last;
	# match already found to this target
      }

      $is_reciprocal++;
    }
  }

  return @hits ? \@hits : undef;
}

sub print_match_stats {
  my ($self) = @_;
  my $stats = $self->match_stats();

  printf STDERR "match types:\n";
  foreach my $k (sort keys %{$stats}) {
    printf STDERR "  %s: %d\n", $k, $stats->{$k};
  }

}

1;

# LRF support:
# _______________      ______________________________      _______________
#                \____/                              \____/
