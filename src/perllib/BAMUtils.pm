package BAMUtils;
# utility functions for BAM files

use strict;
use Exporter;

use Bio::DB::Sam;

use Configurable;
use MiscUtils qw(get_hash_option dump_die);

@BAMUtils::ISA = qw(Configurable Exporter);
@BAMUtils::EXPORT_OK = qw(
			get_bam_rnm
			);

sub get_bam_rnm {
  my ($f_bam) = @_;

  my $bam = Bio::DB::Sam->new(
			      "-bam" => $f_bam,
			      "-expand_flags" => 1,
			     );
  my $rnm = new ReferenceNameMapper();
  foreach my $chrom ($bam->seq_ids()) {
    # standardize reference names
    $rnm->add_name($chrom);
  }
  return $rnm;
}

1;

# LRF support:
# _______________      ______________________________      _______________
#                \____/                              \____/
