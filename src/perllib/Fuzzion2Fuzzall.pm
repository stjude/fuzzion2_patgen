package Fuzzion2Fuzzall;
# perser/helper for fuzzion2 "fuzzall" report
# MNE 2/2024

use strict;
use Exporter;

use Configurable;
use MiscUtils qw(get_hash_option dump_die);

@Fuzzion2Fuzzall::ISA = qw(Configurable Exporter);
@Fuzzion2Fuzzall::EXPORT_OK = qw();

use MethodMaker qw(
		    file
		    df
		    f_pattern
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
  my $file = $self->file || die "-file";

  #
  #  unique-ify header names:
  #

  my $df = new DelimitedFile(
			     "-file" => $file,
			     "-headers" => \&fuzzall_header_fix,
			     );
  # override raw file header with unique-ified headers

  my @hits = grep {/^fuzzall/} @{$df->headers_raw};
  die "can't find fuzzall column" unless @hits == 1;
  $self->f_pattern($hits[0]);
  # column name containing pattern ID (contains fuzzall version #)

  $self->df($df);
}

sub fuzzall_header_fix {
  # STATIC, callback to preprocess raw headers.
  # if this has already been done, no effect.
  my ($h) = @_;
  my %rename = map {$_, 1} qw(
			       min
			       median
			       mean
			       max
			    );

  my $key;
  for (my $i = 0; $i < @{$h}; $i++) {
    if ($rename{$h->[$i]}) {
      $h->[$i] = $key . "_" . $h->[$i];
    } else {
      $key = $h->[$i];
    }
  }
}

sub next {
  my ($self) = @_;
  my $row = $self->df->get_hash();
  if ($row) {
    my @hits = split /,\s*/, $row->{"ID list"};
    my @ids;
    foreach my $hit (@hits) {
      $hit =~ /^(\S+)\((\d+)\/(\d+)\)$/ || die "can't parse pattern ID from $hit";
      my %info;
      $info{id} = $1;
      $info{count_distinct} = $2;
      $info{count_strong_plus} = $3;
      push @ids, \%info;
    }
    $row->{IDs_parsed} = \@ids;
  }
  return $row;
}

1;

# LRF support:
# _______________      ______________________________      _______________
#                \____/                              \____/
