package DWIMColumns;
# wrangle differences in column names

use strict;
use Exporter;

use Configurable;
use MiscUtils qw(get_hash_option dump_die);

@DWIMColumns::ISA = qw(Configurable Exporter);
@DWIMColumns::EXPORT_OK = qw();

use MethodMaker qw(
		    canonical
		    alias2canonical
		    canonical2alias
		    auto_lc
		    canonical2local
		  );

sub new {
  my ($type, %options) = @_;
  my $self = {};
  bless $self, $type;
  $self->canonical({});
  $self->canonical2alias({});
  $self->alias2canonical({});
  # all possible/known aliases

  $self->auto_lc(1);
  $self->configure(%options);
  return $self;
}

sub add_aliases {
  my ($self, $canonical, @aliases) = @_;

  my $alias2canonical = $self->alias2canonical || die;
  my $canonical2alias = $self->canonical2alias || die;
  $self->canonical()->{$canonical} = 1;

  my $auto_lc = $self->auto_lc();

  printf STDERR "add aliases for canonical %s: %s\n", $canonical, join ",", @aliases;

  foreach my $alias ($canonical, @aliases) {
    my @add = $alias;
    if ($auto_lc) {
      my $lc = lc($alias);
      push @add, $lc unless $lc eq $alias;
    }
    foreach my $thing (@add) {
      my $existing = $alias2canonical->{$thing};
      die "collsion for alias $thing" if $existing and $existing ne $canonical;
      $alias2canonical->{$thing} = $canonical;

      $canonical2alias->{$canonical}{$alias} = 1;
    }
  }
}

sub get_data {
  my ($self, $row, $f_wanted) = @_;
  my $f_local = $self->canonical2local->{$f_wanted} || die "can't find local field for $f_wanted";
  return $row->{$f_local};

}

sub init_field_map {
  # for a given field set (i.e. a particular tab-delimited file),
  # find local mappings for desired fields
  my ($self, $set_raw) = @_;
  my $set_local;
  if (ref $set_raw eq "HASH") {
#    $set_local = $set_raw;
    $set_local = { map {$_, 1} keys %{$set_raw} };
  } else {
    die "implement other sets e.g. array";
  }

  my %canonical2local;

#  printf STDERR "local set: %s\n", join ",", sort keys %{$set_local};
  foreach my $f_canonical (keys %{$self->canonical()}) {
    my @aliases = sort keys %{$self->canonical2alias->{$f_canonical}};
#    printf STDERR "aliases for %s: %s\n", $f_canonical, join ",", @aliases;
    foreach my $alias (@aliases) {
      if ($set_local->{$alias}) {
#	printf STDERR "%s => local %s\n", $f_canonical, $alias;
	$canonical2local{$f_canonical} = $alias;
	last;
      }
    }
#    printf STDERR "can't find local entry for %s\n", $f_canonical unless $canonical2local{$f_canonical};

  }
  # TO DO: die if some requested fields not found?

#  dump_die(\%canonical2local, "field map", 1);
  $self->canonical2local(\%canonical2local);
}

sub populate_fields {
  # populate desired field names with data from local mapping
  my ($self, $row) = @_;
  my $canonical2local = $self->canonical2local();
  foreach my $f_canonical (keys %{$canonical2local}) {
    unless (defined $row->{$f_canonical}) {
      # only if column not already populated
      $row->{$f_canonical} = $row->{$canonical2local->{$f_canonical}};
    }
  }
  dump_die($row, "X") unless $row->{sample};
}

1;

# LRF support:
# _______________      ______________________________      _______________
#                \____/                              \____/
