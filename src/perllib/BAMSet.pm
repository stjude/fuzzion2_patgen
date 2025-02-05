package BAMSet;
# index a set of BAM files by SJ sample ID

use strict;
use Exporter;

use File::Basename;

use Configurable;
use FileUtils qw(read_simple_file);
use MiscUtils qw(get_hash_option dump_die);

@BAMSet::ISA = qw(Configurable Exporter);
@BAMSet::EXPORT_OK = qw();

use MethodMaker qw(
    file
    list
    index
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
  # load/index BAM list
  my ($self, %options) = @_;
  my $list;
  if ($list = $self->list()) {
    # pre-parsed arrayref
  } elsif (my $f = $self->file()) {
    $list = read_simple_file($f);
  } else {
    die;
  }

  my %index;
  foreach my $bam (@{$list}) {
    my $bn = basename($bam);
    $bn =~ s/\.bam$//i || die;
    push @{$index{$bn}}, $bam;
#    printf STDERR "save %s %s\n", $bn, $bam;
  }
  $self->index(\%index);
}

sub find {
  my ($self, $sample) = @_;
  return $self->index->{$sample};
}

1;

# LRF support:
# _______________      ______________________________      _______________
#                \____/                              \____/
