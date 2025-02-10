package UniqueAlignmentID;
# describe me

use strict;
use Exporter;

use Configurable;
use MiscUtils qw(get_hash_option dump_die);

@UniqueAlignmentID::ISA = qw(Configurable Exporter);
@UniqueAlignmentID::EXPORT_OK = qw();

use MethodMaker qw(
	unique_counter
        unique_ids
		  );

sub new {
  my ($type, %options) = @_;
  my $self = {};
  bless $self, $type;
  $self->configure(%options);
  $self->unique_counter({});
  $self->unique_ids({});
  return $self;
}

sub get_unique_id {
  my ($self, $a) = @_;
  my $unique_counter = $self->unique_counter();
  my $unique_ids = $self->unique_ids();
  my $strand = $a->strand() == 1 ? "+" : "-";
  my $bn = $a->query->seq_id() . ":" . $strand;
  my $counter = $unique_counter->{$bn}++;
  my $unique_read_id = $bn . ($counter ? $counter : "");
#  printf STDERR "bn:%s final:%s\n", $bn, $unique_read_id;
  die "duplicate" if $unique_ids->{$unique_read_id};
  $unique_ids->{$unique_read_id} = 1;
  return $unique_read_id;
}

1;

# LRF support:
# _______________      ______________________________      _______________
#                \____/                              \____/
