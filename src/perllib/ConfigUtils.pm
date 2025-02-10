package ConfigUtils;

use strict;
use Carp qw(confess);

use Configurable;
use MiscUtils qw(dump_die);
use TdtConfig;
use Exporter;

@ConfigUtils::ISA = qw(Configurable Exporter);
@ConfigUtils::EXPORT_OK = qw(
config_or_manual
);

sub config_or_manual {
  # for a given item, use either manually-specified value, or
  # a value from the config system.
  my %options = @_;
  my $v = $options{"-manual"};
  unless (defined $v) {
    my $config_type = $options{"-config-type"};
    # genome, app
    my $config_name = $options{"-config-name"};

    if ($config_type and $config_name) {
      my $config = $ConfigUtils::CONFIG_CACHE{$config_type}{$config_name};
      unless ($config) {
	$config = TdtConfig::readConfig($config_type, $config_name) || die "can't find config for $config_type $config_name";
	$ConfigUtils::CONFIG_CACHE{$config_type}{$config_name} = $config;
      }
      my $p = $options{"-parameter"} || die "-parameter";
      #    dump_die($config);
      $v = $config->{$p};
    } else {
      confess("config type/name not specified");
    }
  }
  return $v;
}

1;

# LRF support:
# _______________      ______________________________      _______________
#                \____/                              \____/
