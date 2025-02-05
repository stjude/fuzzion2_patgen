package BLASTer;
# BLAST wrapper

use strict;
use Carp qw(confess);

use Bio::SearchIO;
use File::Copy;
use Exporter;

use Digest::MD5;

use Configurable;
use TemporaryFileWrangler;
use WorkingFile;
use FileUtils qw(newer_than);

@BLASTer::ISA = qw(Configurable Exporter);
@BLASTer::EXPORT_OK = qw(
			  dustmasker
		       );

use MethodMaker qw(
tfw
tempfile_base
outfile
null_mode
verbose
strand
ungapped
perc_identity
word_size
dust

output_format
blast_2_sequences

blast_binary
gapopen
gapextend
max_target_seqs
max_hsps

num_alignments
		  );

my @BLAST_INDEX_SUFFIXES = qw(
nhr
nin
nog
nsd
nsi
nsq
);

sub new {
  my ($type, %options) = @_;
  my $self = {};
  bless $self, $type;
  my $tfw = new TemporaryFileWrangler();
#  $tfw->verbose(1);
  $self->tempfile_base($tfw->get_tempfile("-glob" => 1));
  $self->tfw($tfw);
  $self->blast_binary("blastn");
  $self->configure(%options);
  return $self;
}

sub blast {
  my ($self, %options) = @_;

  my $blast2_mode = $self->blast_2_sequences();

  my $output_format = $self->output_format() || "tabular";

  my ($format_blast, $format_parser);

  if ($output_format eq "pairwise") {
    # default
    $format_blast = 0;
    $format_parser = "blast";
    # works?
  } elsif ($output_format eq "tabular") {
    $format_blast = 6;
    $format_parser = "blasttable";
  } elsif ($output_format eq "xml") {
    $format_blast = 5;
    $format_parser = "blastxml";
  } else {
    die "output_format currently limited to tabular/xml";
  }

  my $parser;

  my $verbose = $self->verbose() || $ENV{DEBUG_BLAST} || 0;

  if ($self->null_mode) {
    $parser = Bio::SearchIO->new();
    # not sure this will fly
  } else {
    my $db_fa = $self->fa_setup(($options{"-database"} || die "-database"), "database");
    my $query_fa = $self->fa_setup(($options{"-query"} || die "-query"), "query");
    my $outfile = $self->get_tempfile("-append" => ".blast");
    $self->outfile($outfile);

    printf STDERR "db=%s query=%s out=%s\n", $db_fa, $query_fa, $outfile if $verbose;
    if ($verbose >= 2) {
      foreach my $ref (["query", $query_fa],
		       ["db", $db_fa]) {
	my ($label, $file) = @{$ref};
	printf STDERR "%s FASTA:\n", $label;
	if (-s $file > 5000) {
	    printf STDERR "  (larger file, not printing)\n";
	} else {
	    open(DEBUGTMP, $file) || die;
	    while(<DEBUGTMP>) {
		print STDERR $_;
	    }
	    close DEBUGTMP;
	}
      }

      if ($verbose >= 3) {
	my $f_query_tmp = sprintf 'query_%s.fa', find_fasta_id($query_fa);
	my $f_db_tmp = sprintf 'db_%s.fa', find_fasta_id($db_fa);
	printf STDERR "copied BLAST inputs to %s / %s\n", $f_query_tmp, $f_db_tmp;
	unlink(glob($f_query_tmp . "*"));
	unlink(glob($f_db_tmp . "*"));
	# will unlinks prevent occasional crashes during indexing step?
	copy($query_fa, $f_query_tmp) || die;
	copy($db_fa, $f_db_tmp) || die;
	index_blast_database($f_db_tmp) unless $blast2_mode;
      }
    }

    #
    #  index the database if necessary:
    #
    index_blast_database($db_fa) unless $blast2_mode;

    my @cmd = $self->blast_binary();
    my $protein_mode = $cmd[0] =~ "blastp";

    push @cmd, "-query", $query_fa;
    if ($blast2_mode) {
	push @cmd, "-subject", $db_fa;
    } else {
	push @cmd, "-db", $db_fa;
    }
    push @cmd, "-out", $outfile;
    push @cmd, "-outfmt", $format_blast;
    # TO DO: different formats, maybe w/hash lookup for BioPerl parser below
    push @cmd, "-ungapped" if $self->ungapped();

    push @cmd, "-gapopen " . $self->gapopen() if $self->gapopen();
    push @cmd, "-gapextend " . $self->gapextend() if $self->gapextend();
    push @cmd, "-max_hsps " . $self->max_hsps() if $self->max_hsps();
    push @cmd, "-max_target_seqs " . $self->max_target_seqs() if $self->max_target_seqs();

    unless ($protein_mode) {
      #
      # blastn-specific options:
      #
      if (my $strand = $self->strand) {
	die "unknown strand value $strand, must be plus, minus, or both" unless $strand eq "plus" or $strand eq "minus" or $strand eq "both";
	push @cmd, "-strand", $strand;
      }

      if (my $perc_identity = $self->perc_identity) {
	push @cmd, "-perc_identity", $perc_identity;
      }

      if (my $word_size = $self->word_size) {
	push @cmd, "-word_size", $word_size;
      }

      if (my $dust = $self->dust()) {
	push @cmd, "-dust", $dust;
      }
    }

    my $cmd = join " ", @cmd;

    printf STDERR "running %s\n", $cmd if $verbose;
    if ($verbose >= 3) {
      open(BLASTTMP, ">blast_debug.sh") || die;
      printf BLASTTMP "%s\n", $cmd;
      close BLASTTMP;
    }
    system($cmd);

    die sprintf "error running %s, exit %d", $cmd, $? if $?;
    die "where is outfile $outfile" unless -f $outfile;

    if ($verbose >= 2) {
      printf STDERR "BLAST output:\n";
      open(BLASTTMP, $outfile) || die;
      while (<BLASTTMP>) {
	print STDERR $_;
      }
      close BLASTTMP;
    }

    $parser = Bio::SearchIO->new(-file => $outfile, -format => $format_parser);
  }

  return $parser;
}

sub fa_setup {
  my ($self, $thing, $tag) = @_;
  my $fa_file;
  die unless $tag;
  if (ref $thing) {
    # hash of sequences
    $fa_file = $self->get_tempfile("-append" => sprintf ".%s.fa", $tag);
    write_fasta_set("-file" => $fa_file, "-reads" => $thing);
  } else {
    if (-f $thing) {
      # pre-built filename
      $fa_file = $thing;
    } else {
      # single sequence string
      $fa_file = $self->get_tempfile("-append" => sprintf ".%s.fa", $tag);
      my %reads;
      $reads{query} = $thing;
      write_fasta_set("-file" => $fa_file, "-reads" => \%reads);
    }
  }
  return $fa_file;
}

sub write_fasta_set {
  # hack: should probably figure out how to do w/SeqIO but this is so simple...
  my (%options) = @_;
  my $file = $options{"-file"} || die;
  my $reads = $options{"-reads"} || die;

  if (my @old = glob($file . "*")) {
    unlink @old;
  }

  my $wf = new WorkingFile($file);
  my $fh = $wf->output_filehandle();

  foreach my $id (sort keys %{$reads}) {
    my $sequence = $reads->{$id};
    # QUESTION: mask *?  how does blat deal with these?
    printf $fh ">%s\n%s\n", $id, $sequence;
  }
  $wf->finish();
}

sub get_tempfile {
  my ($self, %options) = @_;
  my $append = $options{"-append"} || die "-append";
  my $base = $self->tempfile_base() || die;
  return $base . $append;
  # glob mode will clean up these varieties
}

sub index_blast_database {
  my ($db_fa) = @_;

  my $needed = 1;

  my @try = (
	     ".nhr",
	     ".00.nhr",
	    );
  # variations on BLAST index file names
  my $idx_file;
  foreach my $suffix (@try) {
    $idx_file = $db_fa . $suffix;
    last if -f $idx_file;
  }

  if (-f $idx_file) {
    $needed = newer_than($db_fa, $idx_file) ? 1 : 0;
    # rebuild if source file is newer than the index file
  } else {
    $needed = 1;
  }

  if ($needed) {
    my $cmd = sprintf 'makeblastdb -in %s -parse_seqids -dbtype nucl > /dev/null', $db_fa;
    system $cmd;
    confess("ERROR: $cmd exited with $?: $!") if $?;
  }
}

sub find_fasta_id {
  # STATIC
  my ($fn) = @_;

  my $md5 = new Digest::MD5();
  open(FAIDTMP, $fn) || die;
  my $first = <FAIDTMP>;
  my $id = "unknown";
  $id = $1 if $first =~ /^>(\S+)/;
  $id =~ s/\W/_/g;
  $md5->add($first);

  while (<FAIDTMP>) {
    $md5->add($_);
  }

  return join "_", $id, $md5->hexdigest();
}

sub dustmasker {
  # STATIC
  # run dustmasker on a single sequence
  my %options = @_;
  my $sequence = $options{"-sequence"} || die;
  my $as_string = $options{"-as-string"};
  my $offset = $options{"-offset"} || 0;

  my $cmd = sprintf 'printf ">query\n%s\n"|dustmasker|', $sequence;
  # hacky but avoids a tempfile

  open(DM, $cmd) || die;
  my @intervals;
  while (<DM>) {
    next if /^>/;
    chomp;
    my @f = split /\s+\-\s+/, $_;
    die unless @f == 2;
    foreach (@f) {
      $_++;
      # convert intervals to 1-based
      $_ += $offset;
    }
    if ($as_string) {
      push @intervals, join "-", @f;
    } else {
      push @intervals, [ @f ];
    }
  }

  return @intervals ? \@intervals : undef;
}

1;

# LRF support:
# _______________      ______________________________      _______________
#                \____/                              \____/
