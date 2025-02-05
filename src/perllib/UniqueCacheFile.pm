package UniqueCacheFile;
# create a cache file name from a set of values uniquely identifying a dataset.
# MNE 9/2014
#
# my $ucf = new UniqueCacheFile("-prefix" => "whatever");
#   # optional prefix, helpful for uniquely identifying these files
# $ucf->add($something1);
# $ucf->add($something2);
# $ucf->add(@list_of_stuff);
# my $outfile = $ucf->get_file();
#
# TO DO:
# - exported single-subroutine version
# - uniqueness checking, if desired

use strict;
use Carp qw(confess);

use Digest::MD5 qw(md5_hex);
use File::Basename;

use Configurable;
use Exporter;

@UniqueCacheFile::ISA = qw(Configurable Exporter);
@UniqueCacheFile::EXPORT_OK = qw();

use MethodMaker qw(
  elements
  prefix
auto_basename
unique
tracker
md5_only
		  );

sub new {
  my ($type, %options) = @_;
  my $self = {};
  bless $self, $type;
  $self->auto_basename(1);
  $self->unique(0);
  $self->tracker({});
  $self->reset();
  $self->configure(%options);
  return $self;
}

sub reset {
  my ($self) = @_;
  $self->elements([]);
}

sub add {
  # add elements required to make this output unique
  my ($self, @things) = @_;
  my $elements = $self->elements() || die;
  foreach my $thing (@things) {
    if (my $ref = ref $thing) {
      if ($ref eq "ARRAY") {
	push @{$elements}, @{$thing};
      } else {
	die "unhandled ref $ref";
      }
    } else {
      push @{$elements}, $thing;
    }
  }
}

sub get_md5 {
  my ($self) = @_;
  return md5_hex(@{$self->elements});
}

sub get_file {
  my ($self, %options) = @_;
  my $retry = $options{"-retry"};

  my @pretty;
  push @pretty, $self->prefix if $self->prefix();

  my $elements = $self->elements();
  die unless @{$elements};

  unless ($self->md5_only) {
    my $auto_basename = $self->auto_basename;
    foreach my $element (@{$elements}) {
      my $s = $element;
      if ($auto_basename and $s =~ /\//) {
	# if an element looks like a filename, only use the basename
	# portion as a printable element.  Prevents very long outfile names.
	# Protected provided against duplicate basenames via MD5.
	$s = basename($s);
      }
      push @pretty, $s;
    }
  }
  push @pretty, $self->get_md5();
  foreach (@pretty) {
    s/\W/_/g;
  }

  my $file = join "_", @pretty;

  if ($self->unique()) {
    # track filenames to ensure files are unique
    # - drawback: can only call this code once!
    confess "duplicate outfile $file" if $self->tracker->{$file} and !$retry;
    # retry mode = we know we're asking again
    $self->tracker->{$file} = 1;
  }

  return $file;
}

1;

# LRF support:
# _______________      ______________________________      _______________
#                \____/                              \____/
