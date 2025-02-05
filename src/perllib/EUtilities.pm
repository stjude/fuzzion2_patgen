package EUtilities;
# "lite" e-utilities implementation
# http://www.ncbi.nlm.nih.gov/books/NBK25499/
# MNE 10/2014

use strict;

use Configurable;
use MiscUtils qw(split_list dump_die shell_cmd);
use UniqueCacheFile;

use LWP::Simple;
use HTTP::Request;
use LWP::UserAgent;
use URI::URL;

use Exporter;

use constant MAINTAINER_EMAIL => 'Michael.Edmonson@stjude.org';
use constant REQUEST_SLEEP_TIME => 1;
# In order not to overload the E-utility servers, NCBI recommends that
# users post no more than three URL requests per second and limit
# large jobs to either weekends or between 9:00 PM and 5:00 AM Eastern
# time during weekdays.

use constant URL_BASE => 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/';

@EUtilities::ISA = qw(Configurable Exporter);
@EUtilities::EXPORT_OK = qw();

use MethodMaker qw(
        database
	retmode
	rettype
	retmax

	cache
        max_query_records
        ucf
	result_files
	retry_count
	all_files
		  );

sub new {
  my ($type, %options) = @_;
  my $self = {};
  bless $self, $type;
  $self->cache(1);
#  $self->max_query_records(100);
  $self->max_query_records(1000);
  $self->retmax(500);
#  $self->max_query_records(20);
  $self->ucf(new UniqueCacheFile("-unique" => 1,
				 "-md5_only" => 1,
				 "-prefix" => "eutilities",
	     ));
  $self->retry_count(5);
  $self->configure(%options);
  return $self;
}

sub fetch_genbank {
  # fetch a set of genbank records by accession
  my ($self, %options) = @_;
  $self->database("nucleotide");
#  $self->database("nuccore");
  $self->retmode("gb");
#  $self->rettype("text");
  $self->rettype("gb");
  return $self->fetch(%options);
}

sub fetch {
  # see http://www.ncbi.nlm.nih.gov/books/NBK25498/#chapter3.Application_3_Retrieving_large

  my ($self, %options) = @_;
  my $ids = $options{"-ids"} || die "-ids";
  # list of query items
  my $retry_count = $self->retry_count();

  if (0) {
    print STDERR "****DEBUG: truncated list!\n";
    $ids = [ @{$ids}[0..10] ];
  }

  my $database = $self->database() || die "database";
  my $retmax = $self->retmax || die;
  my $retmode = $self->retmode || die "retmode";
  my $rettype = $self->rettype || die "rettype";
  my $retmax_query = $self->max_query_records();

  my %pairs_base;
  $pairs_base{db} = $database;
  $pairs_base{tool} = "EUtilities.pm";
  $pairs_base{email} = MAINTAINER_EMAIL;
  $pairs_base{usehistory} = "y";
  # use history server so we don't have to fetch and resubmit IDs

  my $query_url = URL_BASE . "esearch.fcgi";
  my $fetch_url = URL_BASE . "efetch.fcgi";

  my @result_files;
  $self->result_files(\@result_files);

  $self->all_files([]);
  # all cache files, including results as well as Entrez Utilities
  # intermediate files

  foreach my $set (split_list($ids, $retmax_query)) {
    my %pairs_q = %pairs_base;
    $pairs_q{term} = join ",", @{$set};
    $pairs_q{retmax} = $retmax_query;
    # 2/2019: now needed, otherwise default is now apparently 20

    my ($web, $key, $count);
    my $force = 0;

    while (1) {
      my $output_ref = $self->polite_request(
					     "-url" => $query_url,
					     "-pairs" => \%pairs_q,
					     "-content" => 1,
					     "-post" => 1,
					     "-force" => $force,
					    );

      $web = $1 if ($$output_ref =~ /<WebEnv>(\S+)<\/WebEnv>/);
      $key = $1 if ($$output_ref =~ /<QueryKey>(\d+)<\/QueryKey>/);
      $count = $1 if ($$output_ref =~ /<Count>(\d+)<\/Count>/);

      if ($count != scalar @{$set}) {
	# ERROR: sometimes API returns no records, or a count fewer
	# than we asked for
	my $msg = sprintf "EUtilities didn't return requested record count, requested=%d, got=%d", scalar(@{$set}), $count;

	if ($retry_count) {
	  printf STDERR "%s, retrying (%d attempts left)...\n", $msg, $retry_count;
	  $force = 1;
	  $retry_count--;
	} else {
	  # serious problem, better to crash than return data we know is bad
	  die "FATAL ERROR: $msg";
	}
      } else {
	last;
      }
    }

    for (my $retstart = 0; $retstart < $count; $retstart += $retmax) {
      my %pairs_r = %pairs_base;
      $pairs_r{WebEnv} = $web || die;
      $pairs_r{query_key} = $key || die;
      $pairs_r{retstart} = $retstart;
      $pairs_r{retmax} = $retmax;
#      $pairs_r{RetMax} = $retmax;
      $pairs_r{retmode} = $retmode;
      $pairs_r{rettype} = $rettype;

      my $of = $self->polite_request("-url" => $fetch_url,
				     "-pairs" => \%pairs_r
				    );
      push @result_files, $of;
    }
  }

  return \@result_files;
}

sub polite_request {
  # submit query/fetch requests to NCBI servers, sleeping between
  # submissions and and caching results
  my ($self, %options) = @_;
  my $url = $options{"-url"} || die;
  my $pairs = $options{"-pairs"};
  my $ucf = $self->ucf();
  my $force = $options{"-force"};

  $ucf->reset();
  $ucf->add($url);
  my @pairs;
  if ($pairs) {
    @pairs = map {$_, $pairs->{$_}} sort keys %{$pairs};
    $ucf->add(\@pairs)
  }

  my $cache_file = $ucf->get_file("-retry" => $force);
  my $need = not(-s $cache_file);
  $need = 1 if $force;

  if ($need) {
    my $tmp = $cache_file . ".tmp";

    if ($pairs) {
      # pairs of values given: use POST
      # (I think GET will eventually break for very long URLs)
      if ($options{"-post"}) {
	# might need for original query as might be too long for GET
	my $url_o = new URI::URL($url);

	my $request = new HTTP::Request(POST => $url);
	$request->content_type("application/x-www-form-urlencoded");
	$url_o->query_form(@pairs);
	$request->content($url_o->equery());

	my $ua = new LWP::UserAgent();
	my $response = $ua->request($request, $tmp);
	die "error " . $response->error_as_HTML() if $response->is_error;
      } else {
	# larger requests seem to get truncated while fetching results, ???
	# use curl to copy URL instead
	my $url_o = new URI::URL($url);
	$url_o->query_form(@pairs);

	my $tries = 5;

	while ($tries) {
	  my $cmd = sprintf 'curl -o %s "%s"', $tmp, $url_o->as_string();
	  system $cmd;
	  if ($?) {
	    # fail
	    unlink $tmp;
	    if ($tries) {
	      printf STDERR "command %s failed with %s, retrying (%d attempts left)...\n", $cmd, $?, $tries--;
	      sleep 1;
	    } else {
	      die "$cmd failed with $? after retries";
	    }
	  } else {
	    # OK
	    last;
	  }
	}
      }
    } else {
      getstore($url, $tmp);
    }

    die "$tmp missing or empty" unless -s $tmp;
    rename($tmp, $cache_file) || die "can't rename $tmp to $cache_file";
    # only save cache file if download complete

    die "ERROR fetching $url" unless -s $cache_file;
    sleep REQUEST_SLEEP_TIME;
  }

  push @{$self->all_files()}, $cache_file;

  my $result;
  if ($options{"-content"}) {
    local $/ = undef;
    open(EUTMP, $cache_file) || die;
    my $blob = <EUTMP>;
    $result = \$blob;
    close EUTMP;
  } else {
    $result = $cache_file;
  }
  return $result;
}

sub cleanup {
  # remove all cache files:
  #  - intermediate/Entrez Utilities control files
  #  - results files
  my ($self) = @_;
  my $all_files = $self->all_files();
  unlink(@{$all_files}) if @{$all_files};
}

1;

# LRF support:
# _______________      ______________________________      _______________
#                \____/                              \____/
