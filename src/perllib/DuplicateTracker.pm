package DuplicateTracker;
# manage assertions of duplicate relationships
# MNE 7/2020

use strict;
use Exporter;

use Configurable;
use MiscUtils qw(get_hash_option dump_die);

@DuplicateTracker::ISA = qw(Configurable Exporter);
@DuplicateTracker::EXPORT_OK = qw(run_dt_tests);

use MethodMaker qw(
	asserted
	old2new
	verbose
		  );

sub new {
  my ($type, %options) = @_;
  my $self = {};
  bless $self, $type;
  $self->configure(%options);
  $self->reset();
  return $self;
}

sub reset {
  my ($self) = @_;
  $self->asserted({});
  $self->old2new({});
}

sub mark_duplicate {
  my ($self, %options) = @_;
  my $id_from = $options{"-from"} || die "-from";
  my $id_to = $options{"-to"} || die "-to";
  printf STDERR "mark %s dup of %s\n", $id_from, $id_to;
  my $asserted = $self->asserted();
  die "ixnay: duplicate assertion of $id_from" if $asserted->{$id_from};
  $asserted->{$id_from} = 1;
  # problem e.g. if we reject a circular reference and then change the
  # destination:
  #
  # A obsoleted to B
  # B obsoleted to C
  #   - both A and B point to C
  # C obsoleted to A
  #   - detect that this is a circular link and reject: A and B are
  #     suppressed, C survives
  # D obsoleted to A
  #   - D now also points to C
  # A obsoleted to D
  #   - failure?  A now points to D, rather than C, what now is the state of C?
  #   => require that a given key can only be obsoleted once?
  my $old2new = $self->old2new();

  if (my $redirect = $old2new->{$id_to}) {
    # target has already been marked a duplicate
    $id_to = $redirect;
  }
  $old2new->{$id_from} = $id_to;

  while (my ($k, $v) = each %{$old2new}) {
    # update any previous assertions that pointed to $id_from
    # to also now point to $id_to
    $old2new->{$k} = $id_to if $v eq $id_from;
  }
}

sub get_final_id {
  my ($self, $from_id) = @_;
  return $self->old2new->{$from_id} || $from_id;
  # if not marked, return raw ID (i.e. NOP)
}

sub is_duplicate {
  my ($self, $id) = @_;
  my $id_final = $self->get_final_id($id);
  return $id_final ne $id;
}

sub run_dt_tests {
  # STATIC
  my $dt = new DuplicateTracker();
  $dt->verbose(1);

  $dt->reset();
  $dt->mark_duplicate("-from" => "A", "-to" => "B");
  $dt->mark_duplicate("-from" => "B", "-to" => "A");
  die unless $dt->get_final_id("A") eq "B";
  die unless $dt->get_final_id("B") eq "B";
  # both should point to B

  $dt->reset();
  $dt->mark_duplicate("-from" => "A", "-to" => "B");
  $dt->mark_duplicate("-from" => "C", "-to" => "A");
  die unless $dt->get_final_id("A") eq "B";
  die unless $dt->get_final_id("C") eq "B";
  # A and C should both point to B

  $dt->reset();
  $dt->mark_duplicate("-from" => "A", "-to" => "B");
  $dt->mark_duplicate("-from" => "B", "-to" => "C");
  die unless $dt->get_final_id("A") eq "C";
  die unless $dt->get_final_id("B") eq "C";
  # A and B should both point to C

  $dt->reset();
  $dt->mark_duplicate("-from" => "A", "-to" => "B");
  $dt->mark_duplicate("-from" => "B", "-to" => "C");
  $dt->mark_duplicate("-from" => "C", "-to" => "A");
  die unless $dt->get_final_id("A") eq "C";
  die unless $dt->get_final_id("B") eq "C";
  die unless $dt->get_final_id("C") eq "C";
  # C now references itself because A already (indirectly) resolves to C
  # this is essentially a NOP

#  $dt->mark_duplicate("-from" => "A", "-to" => "D");
  # illegal: how to then resolve earlier updates?
  #   A is a dup of B
  #   B is a dup of C
  #
  # maybe a more sophisticated/recursive implementation could deal with this?

}



1;

# LRF support:
# _______________      ______________________________      _______________
#                \____/                              \____/
