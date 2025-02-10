package GenBankCache;
# centralized cache of (versioned) GenBank records.
#
# This is desirable because Entrez Utilities retrieval can be slow and
# glitchy, and the system will also lock out users who overload the
# system with queries.  Multi-threaded operations invoking queries can
# easily trigger a lockout, whereas using cached copies of records
# will not.

# MNE 8/2022
#
# TO DO:
# - option to delete records older than XXX in case updates are available

use strict;
use Exporter;
use File::Path;
use File::Basename;

use Bio::SeqIO;

use FileUtils qw(universal_open);
use Configurable;
use EUtilities;
use MiscUtils qw(get_hash_option dump_die);

@GenBankCache::ISA = qw(Configurable Exporter);
@GenBankCache::EXPORT_OK = qw();

use MethodMaker qw(
	cache_root
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
  unless ($self->cache_root()) {
    my $top_dir = $ENV{HOME} || die "no HOME variable";
    die "where is $top_dir" unless -d $top_dir;
    my $dir = sprintf '%s/.genbank_cache', $top_dir;
    mkpath($dir) unless -d $dir;
    die "where is $dir" unless -d $dir;
    $self->cache_root($dir);
  }
}

sub get {
  #
  # ensure cache is populated for a set of records
  #
  my ($self, $spec, %options) = @_;
  # TO DO: options for catting results into parser stream, etc.
  my $set;
  if (my $ref = ref $spec) {
    if ($ref eq "ARRAY") {
      $set = $spec;
    } else {
      die "unhandled ref type $ref";
    }
  } else {
    # single accession
    $set = [ $spec ];
  }

  my @need;
  foreach my $acc (@{$set}) {
    my $fn = $self->get_cache_file($acc);
    if ($fn and -s $fn) {
      # cache hit, nothing to do
    } else {
      # need to fetch this record
      push @need, $acc;
    }
  }

  if (@need) {
    #
    # if any missing, query Entrez
    #
    my $eu = new EUtilities();
    my $result_files = $eu->fetch_genbank("-ids" => \@need);
    # TO DO: share this code block, similar elsewhere

    foreach my $file (@{$result_files}) {
      my $fh = universal_open($file);

      my @queue;

      my $flush = sub {
	#
	#  write single GenBank record to cache file:
	#
	if (@queue) {
	  my $acc;
	  foreach (@queue) {
	    if (/^VERSION\s+(\S+)/) {
	      $acc = $1;
	      last;
	    }
	  }
	  die "no accession found in file $file: " . join "\n", @queue unless $acc;

	  # write:
	  my $fn = $self->get_cache_file($acc);
	  my $dir = dirname($fn);
	  unless (-d $dir) {
	    mkpath($dir) || die "can't mkpath $dir";
	    die "can't create $dir" unless -d $dir;
	  }
	  open(GBCTMP, ">" . $fn) || die "can't write to $fn: $!";
	  foreach (@queue) {
	    print GBCTMP;
	  }
	  close GBCTMP;
	}
	@queue = ();
      };

      #
      #  extract GenBank records from Entrez Utilities download files:
      #
      while (my $line = <$fh>) {
	&$flush() if $line =~ /^LOCUS/;
	# new record
	push @queue, $line;
      }
      &$flush();
    }

    $eu->cleanup();
  }

  #
  #  return final files:
  #
  my @results;
  foreach my $acc (@{$set}) {
    my $fn = $self->get_cache_file($acc);
    if ($fn and -s $fn) {
      # everything should be available at this point
      push @results, $fn;
    } else {
      die "no local cache file for $acc: $fn";
    }
  }

  return $options{"-single"} ? $results[0] : \@results;
}

sub get_cache_file {
  # local cache filename for a single accession
  my ($self, $acc) = @_;
  my $result;
  my ($record_type, $number, $version);
  if ($acc =~ /^([A-Z]+)_(\d+)\.(\d+)$/) {
    # versioned accession: we know exact name of target file
    ($record_type, $number, $version) = ($1, $2, $3);
  } elsif ($acc =~ /^([A-Z]+)_(\d+)$/) {
    # unversioned accession: might or might not have a cached file available
    ($record_type, $number) = ($1, $2);
  }

  my $dir = sprintf '%s/%s/%s/%s', $self->cache_root(),
    $record_type, substr($number, -2), substr($number, -4);

  if ($version) {
    # versioned accession: we know what the cache filename should be,
    # even if we don't yet have it
    $result = sprintf '%s/%s', $dir, $acc;
  } else {
    # unversioned accession: only return filename if it exists
    my $glob = sprintf '%s/%s*', $dir, $acc;
    my @files = glob($glob);
    if (@files == 1) {
      $result = $files[0];
    } elsif (@files > 1) {
      # if multiple records, return the latest version
      my %suffix2file = map {/\.(\d+)$/ => $_} @files;
      my @ordered = sort {$b <=> $a} keys %suffix2file;
      $result = $suffix2file{$ordered[0]};
    }
  }

  return $result;
}

1;

# LRF support:
# _______________      ______________________________      _______________
#                \____/                              \____/
