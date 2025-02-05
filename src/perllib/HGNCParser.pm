package HGNCParser;
# parser/lookup for HGNC (HUGO) gene symbols & synonyoms

use strict;
use Carp qw(confess);
use Configurable;
use Exporter;

use DelimitedFile;
use MiscUtils qw(dump_die unquote);

@HGNCParser::ISA = qw(Configurable Exporter);
@HGNCParser::EXPORT_OK = qw();

use MethodMaker qw(
	file

	ids

	index_approved
	index_refseq
	index_synonym
	index_prev
	index_all_sym

lite
blacklist_prev
blacklist_synonym
		  );

sub new {
  my ($type, %options) = @_;
  my $self = {};
  bless $self, $type;
  $self->blacklist_prev({});
  $self->blacklist_synonym({});
  $self->configure(%options);
  $self->setup();
  return $self;
}

sub setup {
  my ($self, %options) = @_;
  my $fn = $self->file || die;
  my $df = new DelimitedFile(
			     "-file" => $fn,
			     "-headers" => 1,
			     );

  my %ids;

  my %idx_approved;
  my %idx_nm;
  my %idx_synonym;
  my %idx_prev;
  my %idx_all_sym;

  my %allsym2id;
  # all symbols

  my $lite_mode = $self->lite();
  my @f_lite_wanted = (
		       "HGNC ID",
		       "Approved Symbol",
		       "Previous Symbols",
		       "Synonyms",
		       "RefSeq IDs",
		      );

  while (my $row = $df->get_hash()) {

    if (exists $row->{hgnc_id}) {
      # newer file format.  Convert basic fields to old format
      # for backwards compatibility.

      $row->{"HGNC ID"} = $row->{hgnc_id};
      $row->{"Approved Symbol"} = $row->{symbol};

      $row->{"Previous Symbols"} = convert_new_list_to_old($row->{prev_symbol});
      $row->{"Synonyms"} = convert_new_list_to_old($row->{alias_symbol});
      $row->{"RefSeq IDs"} = convert_new_list_to_old($row->{refseq_accession});
    }

    if ($lite_mode) {
      # to conserve RAM, only keep minimum fields
      my %r;
      @r{@f_lite_wanted} = @{$row}{@f_lite_wanted};
      $row = \%r;
    }

    my $id = $row->{"HGNC ID"} || die;
    die "duplicate $id" if $ids{$id};
    $ids{$id} = $row;
    # store rows by unique ID, indexes point to IDs

    my $approved_sym = $row->{"Approved Symbol"} || die;
    $idx_approved{$approved_sym}{$id} = 1;
    $idx_approved{uc($approved_sym)}{$id} = 1;
    # use same structure as other indexes

    my @all_syms = $approved_sym;

    if (my $refs = $row->{"RefSeq IDs"}) {
      my @rs = split /,\s*/, $refs;
      foreach my $rs (@rs) {
	$idx_nm{$rs}{$id} = 1;
      }
      $row->{refseqs_hash} = { map {$_, 1} @rs };
    }

    foreach my $ref (
		     [ "Previous Symbols", \%idx_prev ],
		     [ "Synonyms", \%idx_synonym ],
		    ) {
      my ($field, $hash) = @{$ref};
      die unless exists $row->{$field};
      if (my $things = $row->{$field}) {
	my @things = split /,\s*/, $things;
	push @all_syms, @things;
	foreach my $thing (@things) {
	  $hash->{uc($thing)}{$id} = 1;
	}
      }
    }
#    printf STDERR "all for %s => %s\n", $approved_sym, join ",", @all_syms;

    foreach my $sym (@all_syms) {
      $idx_all_sym{uc($sym)}{$id} = 1;
    }
  }

  $self->ids(\%ids);
  $self->index_approved(\%idx_approved);
  $self->index_refseq(\%idx_nm);
  $self->index_prev(\%idx_prev);
  $self->index_synonym(\%idx_synonym);
  $self->index_all_sym(\%idx_all_sym);
}

sub find {
  my ($self, %options) = @_;
  my $symbol = $options{"-symbol"};
  my @ids;
  if (defined $symbol) {
    $symbol = uc($symbol);
    my $hit;
    if ($options{"-previous"}) {
      $hit = $self->find_by_index($self->index_prev(), $symbol);
    } elsif ($options{"-approved"}) {
      $hit = $self->find_by_index($self->index_approved(), $symbol);
    } elsif ($options{"-synonym"}) {
      $hit = $self->find_by_index($self->index_synonym(), $symbol);
    } else {
      # any symbol
      $hit = $self->find_by_index($self->index_all_sym(), $symbol);
    }
    die if $hit and not(ref $hit);
    @ids = keys %{$hit} if $hit;
  } else {
    # other types?
    die "-symbol";
  }

  my @results;
  if (@ids) {
    @results = map {$self->ids->{$_} || die "no ID for $_"} @ids;
  }

  return @results ? \@results : undef;
}

sub is_approved {
  my ($self, $sym) = @_;
#  my $status = $self->index_approved->{$sym};
  return $self->find_by_index($self->index_approved, $sym);
}

sub is_previous {
  my ($self, $sym) = @_;
  return $self->find_by_index($self->index_prev, $sym);
}

sub is_synonym {
  my ($self, $sym) = @_;
  return $self->find_by_index($self->index_synonym, $sym);
}

sub find_by_index {
  my ($self, $index, $sym) = @_;
  return $index->{uc($sym)};
  # symbols are always uppercased in index
}

sub prune_synonyms_old {
  my ($self) = @_;

  my $idx_synonym = $self->index_synonym();
  my $idx_approved = $self->index_approved();
  foreach my $sym (sort keys %{$idx_synonym}) {
    my $bad;
    if ($idx_approved->{$sym}) {
      # if a synonym is itself an approved symbol, remove link to other gene.
      # e.g. FAH is both an approved symbol as well as a synonym for FANCA.
#      printf STDERR "%s is both synonym and approved\n", $sym;
      $bad = 1;
    }

    if (keys %{$idx_synonym->{$sym}} > 1) {
      # e.g. 2F1 is ambiguous: maps to both SLC25A5 and KLRG1
#      printf STDERR "ambiguous synonym %s => %s\n", $sym, join ",", sort keys %{$idx_synonym->{$sym}};
      $bad = 1;
    }
    
    if ($bad) {
      delete $idx_synonym->{$sym};
    }
  }
}

sub prune_synonyms {
  my ($self) = @_;
  $self->prune_index("-index" => $self->index_synonym,
		     "-label" => "synonyms",
		     "-blacklist" => $self->blacklist_synonym,
		     "-synonym" => 1
      );
}

sub prune_previous {
  my ($self) = @_;
  $self->prune_index("-index" => $self->index_prev,
		     "-label" => "previous",
		     "-blacklist" => $self->blacklist_prev
      );
  # e.g. TCF3 and TCF7L1: both are approved symbols, but TCF3 is also
  # a previous symbol for TCF7L1
}

sub prune_index {
  my ($self, %options) = @_;
  my $idx = $options{"-index"} || die;
  my $label = $options{"-label"} || die;
  my $blacklist = $options{"-blacklist"} || die;
  my $synonym_mode = $options{"-synonym"};
      
  my $idx_approved = $self->index_approved();
  my $verbose = $ENV{HGNC_VERBOSE};
  foreach my $sym (sort keys %{$idx}) {
    my $bad;
    if ($idx_approved->{$sym}) {
      # if symbol is also an approved symbol, remove link to other gene.
      # e.g. FAH is both an approved symbol as well as a synonym for FANCA.
      printf STDERR "prune %s: %s is also an approved symbol\n", $label, $sym if $verbose;
      $bad = "also_approved";
    }

    if (keys %{$idx->{$sym}} > 1) {
      # e.g. 2F1 is ambiguous: maps to both SLC25A5 and KLRG1
      printf STDERR "prune %s: %s is ambiguous (%s)\n", $label, $sym, join ",", sort keys %{$idx->{$sym}} if $verbose;
      $bad = "ambiguous";
    }

    if ($synonym_mode and $self->index_prev()->{$sym}) {
      # additional scrutiny for synonyms: does a synonym also appear
      # in the previous symbol index?  e.g. MLL2 is a synonym for KMT2B,
      # however it's a PREVIOUS symbol for KMT2D.  In this case prune
      # the synonym entry.
      $bad = "also_previous";
    }
    
    if ($bad) {
      delete $idx->{$sym};
      $blacklist->{$sym} = $bad;
    }
  }
}

sub is_blacklisted_previous {
  my ($self, $sym) = @_;
  return $self->blacklist_prev->{$sym};
}

sub is_blacklisted_synonym {
  my ($self, $sym) = @_;
  return $self->blacklist_synonym->{$sym};
}

sub find_hgnc_record {
  my ($self, %options) = @_;
  my $sym = $options{"-symbol"};
  confess("-symbol") unless defined $sym;
  my $result;
  if (my $thing = $self->is_approved($sym)) {
    # FIX ME: what about alt types??
    my ($hgnc_id) = (keys %{$self->index_approved->{uc($sym)}});;
    # FIX ME: some kind of accessor method rather than requiring manual uc()
    die "approved symbol $sym but no HGNC ID, ?" unless $hgnc_id;
    $result = $self->ids()->{$hgnc_id} || die;
  }
  return $result;
}

sub convert_new_list_to_old {
  # STATIC
  my ($old) = @_;
  return join ", ", split /\|/, unquote($old);
}

sub is_blacklisted {
  my ($self, $sym) = @_;

  my @reasons;
  my $b_prev = $self->is_blacklisted_previous($sym);
  my $b_syn = $self->is_blacklisted_synonym($sym);
  push @reasons, sprintf 'previous=%s', $b_prev if $b_prev;
  push @reasons, sprintf 'synonym=%s', $b_syn if $b_syn;
  return join ",", @reasons;
}

sub get_approved_symbols {
  my ($self) = @_;
  return [ keys %{$self->index_approved()} ];
}

1;

# LRF support:
# _______________      ______________________________      _______________
#                \____/                              \____/
