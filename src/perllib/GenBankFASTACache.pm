package GenBankFASTACache;
# cache FASTA sequences retrieved from GenBank
# MNE 8/2021

use strict;
use Exporter;

use Configurable;
use MiscUtils qw(get_hash_option dump_die);
use FAI qw(build_fai);
use EUtilities;
use GeneSymbolMapper qw(new_gsm_lite);

@GenBankFASTACache::ISA = qw(Configurable Exporter);
@GenBankFASTACache::EXPORT_OK = qw();

use MethodMaker qw(
		    file
		    hgnc

		    gene2acc
		    gsm
		    fai
		    acc_versioned
		    acc_unversioned
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
  my ($self) = @_;
  my $fn = $self->file() || die;

  my $gsm = new_gsm_lite();

  my $fai;
  my %gene2acc;
  my %acc_versioned;
  my %acc_unversioned;

  if (-f $fn) {
    # OK if doesn't exist yet
    open(NCTMP, $fn) || die "can't open $fn: $!";
    while (<NCTMP>) {
      if (/^>/) {
	/^>(\S+)/;
	my $acc = $1;
	/gene=(\S+)/ || die;
	# ugly; maybe store this metadata in different file??
	my $gene = $1;
	push @{$gene2acc{$gene}}, $acc;
	$gsm->add_gene("-gene" => $gene);

	$acc_versioned{$acc} = 1;
	$acc =~ s/\.\d+$//;
	$acc_unversioned{$acc} = 1;
      }
    }
    close NCTMP;

    $fai = new FAI("-fasta" => $fn);
  }

  $self->gene2acc(\%gene2acc);
  $self->acc_versioned(\%acc_versioned);
  $self->acc_unversioned(\%acc_unversioned);
  $self->gsm($gsm);
  $self->fai($fai);
}


sub add_gene {
  # add a sequence for a gene, using HGNC for the lookup sequence.
  # TO DO: use SJPI as well?
  # HGNC is useful because this code is typically used for non-refseq sequences

  my ($self, $gene) = @_;

  die sprintf "ERROR: record aready present for %s in %s\n", $gene, $self->file() if $self->gsm->find($gene);

  my $f_hgnc = $self->hgnc() || die "-hgnc";
  printf STDERR "HGNC: %s\n", $f_hgnc;

  my $df = new DelimitedFile("-file" => $f_hgnc,
			     "-headers" => 1,
			     );
  my $gsm2 = new_gsm_lite();
  my %sym2acc;

  while (my $row = $df->get_hash()) {
    my $sym = $row->{symbol} || die;
    my $acc = $row->{refseq_accession};

    if ($acc) {
      # not always present
      $sym2acc{$sym} = $acc;
      $gsm2->add_gene("-gene" => $sym);
    }
  }

  my $acc = $sym2acc{$gene};
  unless ($acc) {
    if (my $alt_sym = $gsm2->find($gene)) {
      $gene = $alt_sym;
      # HGNC file uses a different primary symbol, use that instead (?)
      $acc = $sym2acc{$gene};
    }
  }

  die "can't find accession for $gene" unless $acc;


  $self->add_accession($acc, $gene);
}

sub add_accession {
  my ($self, $acc, $gene) = @_;
  die unless $acc and $gene;

  if ($self->have_accession($acc)) {
    printf STDERR "already have a record for %s\n", $acc;
    return;
  }

  #
  # fetch FASTA sequence from GenBank:
  #
  my $eu = new EUtilities();
  $eu->database("nucleotide");
  $eu->rettype("fasta");
  $eu->retmode("text");

  my $hits = $eu->fetch("-ids" => [ $acc ]);
  die "error, multiple files" unless @{$hits} == 1;
  my $f_gb_raw = $hits->[0];

  #
  # merge with existing FASTA file, adding gene annotation
  #
  my $f_fasta = $self->file();
  open(IN, $f_gb_raw) || die;
  open(OUT, ">>" . $f_fasta) || die;
  while (<IN>) {
    if (/^>/) {
      chomp;
      printf OUT "%s /gene=%s\n", $_, $gene;
    } else {
      print OUT if /\w/;
      # ignore blank lines (trailing)
    }
  }
  close IN;
  close OUT;

  build_fai($f_fasta);
  # index fasta file

  $self->setup();
  # rebuild internal index

  # TO DO: clean up eutilities cache files?
}

sub has_gene {
   my ($self, $gene) = @_;
   return $self->gsm->find($gene);
}

sub get_sequences_for_gene {
   my ($self, $gene_raw, %options) = @_;
   my @results;
   # TO DO: option to return as hash
   if (my $gene = $self->gsm->find($gene_raw)) {
     my $set = $self->gene2acc->{$gene} || die;
     my $fai = $self->fai;

     foreach my $acc (@{$set}) {
       my $seq_ref = $fai->get_sequence("-id" => $acc);
       my %row;
       $row{accession} = $acc;
       $row{seq_ref} = $seq_ref;
       push @results, \%row;
     }
   }
   return \@results;
}

sub have_accession {
  my ($self, $acc) = @_;
  my $hash_check = $acc =~ /\.\d+$/ ? $self->acc_versioned : $self->acc_unversioned;
  # if request is unversioned, check unversioned. If request is
  # versioned, check versioned.  This lets us support multiple
  # versions if that level of granualarity is used.
  return $hash_check->{$acc};
}

1;

# LRF support:
# _______________      ______________________________      _______________
#                \____/                              \____/
