package BLATClient;
# wrapper to SJ BLAT server
# MNE 8/2013

use strict;
use Carp qw(confess);

use Configurable;
use Exporter;
use FileUtils qw(find_binary newer_than);
use MiscUtils qw(dump_die);
#use misc;
use TdtConfig;

use Bio::SearchIO;

@BLATClient::ISA = qw(Configurable Exporter);
@BLATClient::EXPORT_OK = qw();

use MethodMaker qw(
program
host
port
genome
root
outfile
outfile_literal
retries

format
config_genome
headered

minScore

		  );

my %GENOME2PORT = (
  "hg18" => 50000,
  "hg19" => 50037,
  "mm9" => 50001,
    );
# servers managed by St. Jude / HPCF

sub new {
  my ($type, %options) = @_;
  my $self = {};
  bless $self, $type;
  $self->program("gfClient");
  #James Mc
#2 blat servers have been set up at "clinicalblat:50037" against 
#/clingen/dev/tartan/index/reference/Homo_sapiens/GRCh37-lite/2BIT/GRCh37-lite.2bit 
#for the clinical cluster.
#$self->host("clinicalblat");
  $self->host("sjblat");
  # Mi Zhou:
  # "sjblat" is the correct domain name to use.
  # It resolves to 4 blat servers round-robin.
#  $self->host("sjblat.stjude.org");
  $self->genome("hg19");
  $self->root("/");
  $self->retries(10);
  $self->format("psl");
  # output format
  $self->headered(1);
  $self->configure(%options);
  return $self;
}

sub query {
  my ($self, %options) = @_;
  my $fa = $options{"-fasta"} || die "-fasta";

  my $format = $self->format() || die;

  my $outfile;
  if ($self->outfile_literal) {
    $outfile = $self->outfile() || die;
  } else {
    $outfile = ($self->outfile() || $options{"-outfile"} || $fa) . "." . $format;
  }

  my $options = $options{"-options"} || "";
  $options .= " -nohead" unless $self->headered();
  $options .= sprintf " -minScore=%d", $self->minScore() if defined $self->minScore();

  my $port;
  if (my $cg = $self->config_genome) {
    # new method (preferred)
    my $config_genome = TdtConfig::readConfig('genome', $cg) || die "can't find config for $cg";
    $port = $config_genome->{BLAT_PORT} || die;
    $self->host($config_genome->{BLAT_HOST} || die);
  } elsif (my $genome = $self->genome()) {
    # DEPRECATED
    $port = $GENOME2PORT{$genome} || die "unknown port for $genome";
  } elsif ($self->host and $self->port) {
    # externally specified
    $port = $self->port();
  } else {
    die "specify -genome [deprecated] or config_genome";
  }

  my $program = $self->program();
  find_binary($program, "-die" => 1);
  # die unless required module/app profile is loaded

  my $fmt = sprintf 'out=%s', $format;

  my $cmd = join " ",
  $program,
  $self->host(),
  $port,
  $self->root(),
  $fa,
  $outfile,
  $fmt,
  $options;

  printf STDERR "blat client: %s\n", $cmd if $ENV{VERBOSE};

  my $run_needed = not(-f $outfile) || newer_than($fa, $outfile);

  if ($run_needed) {
    my $max_tries = $self->retries();
    my $try;
    for ($try = 1; $try <= $max_tries; $try++) {
      system $cmd;
      if ($?) {
	if ($try == $max_tries) {
	  confess "ERROR calling $cmd after retries; gfClient not on path or BLAT server down? exit code=$?";
	} else {
#	  my $sleep = $try * 2;
	  my $sleep = $try ** 3;
	  # 1st attempt: 1 second
	  # 2nd attempt: 8 seconds
	  # 3rd attempt: 27 seconds
	  # ...
	  # 10th attempt: ~16 m
	  printf STDERR "error invoking BLAT server (code %s: %s), sleeping %d and retrying...\n", $?, $!, $sleep;
	  sleep $sleep;
	}
      } else {
	# OK
	last;
      }
    }
  }

  my $parser = Bio::SearchIO->new(-file => $outfile, -format => 'psl');
  # ugh: no specific pslx parser?
  return $parser;
}

1;

# LRF support:
# _______________      ______________________________      _______________
#                \____/                              \____/               
