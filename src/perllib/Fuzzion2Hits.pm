package Fuzzion2Hits;
# parser/interator for Fuzzion2 .hits files
# MNE 6/2025

use strict;
use FileHandle;
use Exporter;

use Configurable;
use MiscUtils qw(get_hash_option dump_die);
use WorkingFile;

@Fuzzion2Hits::ISA = qw(Configurable Exporter);
@Fuzzion2Hits::EXPORT_OK = qw();

use MethodMaker qw(
    hits

    header
    header_line

    fh_hits

    set_lines
    row_info

    out
    out_wf

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
  my $f_hits = $self->hits || die "-hits";
  my $fh = new FileHandle();
  $fh->open($f_hits) || die;
#  print STDERR "hits $f_hits\n";
  my $h = <$fh>;
  chomp $h;
  $self->header_line($h);
  my @h = split /\t/, $h;
  $self->header(\@h);
  $self->fh_hits($fh);

  if (my $f_out = $self->out) {
    my $wf = new WorkingFile($f_out);
    $self->out_wf($wf);
    my $fh_out = $wf->output_filehandle;
    printf $fh_out "%s\n", $h;
  }
}

sub next_row {
  my ($self) = @_;
  my $result;

  my $fh = $self->fh_hits || die;
  if (my $l_p = <$fh>) {
    if ($l_p =~ /^read\-pairs /) {
      # "read-pairs": final summary info.  TO DO: parse if wanted?
      if (my $wf = $self->out_wf()) {
	my $fh_out = $wf->output_filehandle();
	print $fh_out $l_p;
      }
    } else {
      die $l_p unless $l_p =~ /^pattern /;

      my $l_r1 = <$fh> || die;
      die unless $l_r1 =~ /^read /;
      my $l_r2 = <$fh> || die;
      die unless $l_r2 =~ /^read /;

      foreach ($l_p, $l_r1, $l_r2) {
	chomp;
      }

      $self->set_lines([$l_p, $l_r1, $l_r2]);

      my @f = split /\t/, $l_p, -1;
      # -1: if last column value blank
      my $h = $self->header();
      die sprintf "column count mismatch: headers=%d row=%d, %s",
	  scalar(@{$h}), scalar @f, $l_p unless @f == @{$h};

      my %r;
      @r{@{$h}} = @f;

      $self->row_info(\%r);

      $result = \%r;
    }
  }

  return $result;
}

sub write_row {
  my ($self) = @_;
  if (my $wf = $self->out_wf()) {
    my $fh_out = $wf->output_filehandle();
    my $list = $self->set_lines();
    foreach my $line (@{$list}) {
      printf $fh_out "%s\n", $line;
    }
  } else {
    die "can't write, not initialized";
  }

}

sub finish {
  my ($self) = @_;
  $self->out_wf()->finish() if $self->out_wf();
}

sub get_pid {
  my ($self) = @_;
  my $pid = $self->row_info()->{$self->header->[0]};
  $pid =~ s/^pattern // || die;
  return $pid;
}

1;

# LRF support:
# _______________      ______________________________      _______________
#                \____/                              \____/
