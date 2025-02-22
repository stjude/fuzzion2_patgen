package DBTools;

@DBTools::ISA = qw(Configurable Exporter);
@DBTools::EXPORT_OK = qw(
  fetch_one
  selectall_hashref
  selectall_colhash
  hash_insert
  hash_update

  get_dbi_lpg
  get_dbi_hg19
  get_dbi_hg18
  get_dbi_hg
  get_dbi_hg_readonly
  get_dbi_cgwb_readonly
  get_dbi_heatmap_readonly

  get_dbi_mm9

  get_dbi_snp
  get_dbi_tcga_brca_indel

get_dbi_gedi
get_dbi_hgmd

  export_query_to_flatfile
parse_pg_pass
parse_pg_service
  );

use strict;

use DBI;
#use Date::Manip;
use Sys::Hostname;
use Carp qw(confess);

use Configurable;
use MiscUtils qw(dump_die);
use TdtConfig;

use MethodMaker qw(
		   dbi

                   is_mysql
                   is_oracle
		  );

sub new {
  my ($type, %options) = @_;
  my $self = {};
  bless $self, $type;
  $self->configure(%options);
  $self->setup();
  return $self;
}

sub grouper {
  # group by specified field and call back
  my ($self, $cmd, $callback, %options) = @_;
  my $dbi = $self->dbi() || die "need dbi";
  my $column = $options{"-column"} || 0;
  
  my $sth = $dbi->prepare($cmd) || die "prep failed";
  $sth->{LongReadLen} = 10000;
  # try to stave off ORA-01406 truncation errors for "long" fields
  $sth->execute() || die "exec failed for $cmd";
  my $last = "";
  my @rows;
  my @row;
  while (@row = $sth->fetchrow_array()) {
 #   printf "row: %s\n", join ", ", @row;
    if ($row[$column] ne $last) {
      &$callback(\@rows) if @rows;
      @rows = ();
    }
    $last = $row[$column];
    push @rows, [ @row ];
  }
  $sth->finish();
  &$callback(\@rows) if @rows;
}

sub fetch_one {
  my ($sth, @args) = @_;
  $sth->execute(@args);
  my $ref = $sth->fetchrow_arrayref();
  if ($ref) {
    while ($sth->fetchrow_arrayref) {};
  }
  return $ref;
}

sub selectall_hashref {
  # return arrayref of hashrefs for each row
  my ($dbi, $cmd, %options) = @_;
  my $sth = $dbi->prepare($cmd) || die;
  $sth->execute() || die "execute of $cmd failed";
  my @results;
  while (my $ref = $sth->fetchrow_hashref) {
    push @results, $ref;
  }
  $sth->finish;
  if ($options{"-names"}) {
    # return ordered column names as well
    return (\@results, $sth->{NAME});
  } else {
    return \@results;
  }
}

sub selectall_colhash {
  # return a single-column query as a hashref
  my ($dbi, $cmd) = @_;
  my $sth = $dbi->prepare($cmd) || die;
  $sth->execute() || die;
  my %results;
  my $ref;
  while ($ref = $sth->fetchrow_arrayref) {
    $results{$ref->[0]} = 1;
  }
  $sth->finish;
  return \%results;
}

sub setup () {
  my ($self) = @_;
  my $dbi = $self->dbi || die "need -dbi";

  my $type = $dbi->get_info(17);
  $self->is_mysql(0);
  $self->is_oracle(0);

  if ($type eq "Oracle") {
    $self->is_oracle(1);
  } elsif ($type eq "MySQL") {
    $self->is_mysql(1);
  } else {
    die "don't know how to initialize for this DBI!";
  }
}

sub standardize_date {
  my ($self, $date_string) = @_;

  if ($self->is_mysql) {
    # With Oracle we can use NLS_DATE_FORMAT and ALTER SESSION to
    # specify date formats during insertion, but MySQL wants
    # dates in ISO format (YYYY-MM-DD).
    if (my $dm = ParseDate($date_string)) {
      # date format parsable
      $date_string = UnixDate($dm, "%Y/%m/%d");
      # reformat
    }
  }
  
  return $date_string;
}

sub dbi_verify {
  my ($var) = @_;
  return ($var and $var->ping() ? 1 : 0);
}

sub get_mysql_attr {
  my (%options) = @_;
  my %attr;
  return \%attr;
}

sub get_dbi_hg_readonly {
  confess("unavailable");
}

sub get_dbi_hg19 {
  confess("unavailable");
}

sub get_dbi_hgmd {
  confess("unavailable");
}

sub get_dbi_hg18 {
  confess("unavailable");
}

sub get_dbi_mm9 {
  confess("unavailable");
}

sub get_dbi_hg {
  confess("unavailable");
}

sub get_dbi_snp {
  confess("unavailable");
}

sub get_dbi_tcga_brca_indel {
  confess("unavailable");
}

sub get_dbi_lpg {
  confess("unavailable");
}

sub get_dbi_cgwb_readonly {
  confess("unavailable");
}


sub get_dbi_heatmap_readonly {
  confess("unavailable");
}

sub hash_insert {
  my (%options) = @_;
  my $dbi = $options{"-dbi"} || die "-dbi";
  my $table = $options{"-table"} || die "-table";
  my $row = $options{"-row"} || die "-row";

  if (my $adf = $options{"-auto-delete"}) {
    my $cmd = sprintf "delete from %s where %s = ?", $table, $adf;
    my $sth_delete = $dbi->prepare($cmd) || die "prep failed for $cmd";
    $sth_delete->execute($row->{$adf} || die "no autodel value for $adf") || die || die "autodelete failed";
  }
  
  my @fields = keys %{$row};
  # sort order doesn't matter but must be consistent
  my $cmd = sprintf 'INSERT INTO %s (%s) VALUES (%s)', $table, 
      join(",", @fields), join (",", split //, "?" x scalar @fields);
  my $sth = $dbi->prepare($cmd) || die "prepare failed for $cmd";
  $sth->execute(@{$row}{@fields}) || die "execute failed!";
}

sub hash_update {
  my (%options) = @_;
  my $dbi = $options{"-dbi"} || die "-dbi";
  my $table = $options{"-table"} || die "-table";
  my $row = $options{"-row"} || die "-row";
  my $index_col = $options{"-index-column"} || die "-index-column";
  my $existing = $options{"-existing"};
  my $not_zero = $options{"-no-zero-update"};
  my $not_empty = $options{"-no-empty-update"};
  my $verbose = $options{"-verbose"};
  $verbose = 1 unless defined $verbose;
  unless ($existing) {
    #
    # get current row from database
    #
    my $sth = $dbi->prepare(sprintf "select * from %s where %s = ?", $table, $index_col) || die;
    $sth->execute($row->{$index_col}) || die;
    my @hits;
    while (my $ref = $sth->fetchrow_hashref()) {
      push @hits, $ref;
    }
    die "record not found" unless @hits;
    die "record ambiguous" if @hits > 1;
    $existing = $hits[0];
  }
  my %allow_zero;
  if (my $az = $options{"-allow-zero"}) {
    %allow_zero = map {$_, 1} @{$az};
  }

  my %different;

  foreach my $col (keys %{$row}) {
#    printf STDERR "%s: %s %s\n", $col, $existing->{$col}, $row->{$col};
    if (($existing->{$col} || "") ne ($row->{$col} || "")) {
      # difference between existing and new values
#      printf STDERR "hash_update(): difference in %s old=%s new=%s\n", $col, ($existing->{$col} || ""), ($row->{$col} || "");
      my $new_value = defined $row->{$col} ? $row->{$col} : "";
      next if $not_empty and $new_value eq "";
#      next if $not_zero and $new_value eq "0";

      if ($new_value eq "0" and $not_zero) {
	next unless $allow_zero{$col};
      }
      
      printf STDERR "hash_update(): difference in %s old=%s new=%s\n",
	$col,
	  (defined $existing->{$col} ? $existing->{$col} : ""),
	    $new_value if $verbose;
#	    (defined $row->{$col} ? $row->{col} : "");
      
      my $cmd = sprintf 'update %s set %s=? where %s=?', $table, $col, $index_col;
      my $sth = $dbi->prepare($cmd) || die;
      $sth->execute($row->{$col}, $row->{$index_col}) || die "update failed";

      if (0 and $col eq "ProjectID") {
	printf STDERR "DEBUG: NOT changing ProjectID from %s => %s\n", $existing->{$col}, $row->{$col};
	next;
      }

      if (0) {
	unless ($col eq "mtime" or $col eq "SampleID" or $col eq "ProjectID" or $col eq "ReadLength" or $col eq "bam_run_date") {
	  die join ",", "CHECK ME:", $cmd, $row->{$col}, $row->{$index_col};
	}
      }

    }
  }

}

sub get_dbi_gedi {
  # http://hc-wiki.stjude.org/display/compbio/Data+Manager+New+Hire+Checklist
  my (%options) = @_;
  my $type = $options{"-type"} || die "-type [research|clinical]";
  my $db_user = getpwuid($<) || die "can't identify login";

  my $db_password = $options{"-password"};
  unless ($db_password) {
    if ($type eq "research") {
      $db_password = $ENV{GEDI_PASSWORD_RESEARCH};
    } elsif ($type eq "clinical") {
      $db_password = $ENV{GEDI_PASSWORD_CLINICAL};
    } else {
      die;
    }
  }

  unless ($db_password) {
    printf STDERR "GeDI password: ";
    $db_password = <STDIN>;
    chomp $db_password;
  }

  my $config = TdtConfig::readConfig("dbconn", "raptr") || die "can't find config for dbconn/raptr";
  my $dsn = $config->{DSN} || die "no DSN";

  my $dbi = DBI->connect($dsn, $db_user,$db_password, {RaiseError => 1, AutoCommit => 0});

  return $dbi;
}

sub export_query_to_flatfile {
  my (%options) = @_;
  my $dbi = $options{"-dbi"} || die "-dbi";
  my $sql;
  my $outfile = $options{"-outfile"};
  my $trim_whitespace = $options{"-trim-whitespace"};
  my $header = $options{"-header"};
  $header = 1 unless defined $header;

  if (my $table = $options{"-table"}) {
    # table or view
    $sql = sprintf 'select * from %s', $table;
    $outfile = sprintf 'export_%s.tab', $table unless $outfile;
  } elsif ($sql = $options{"-sql"}) {
    # ok
  } elsif (my $f_sql = $options{"-sql-file"}) {
    # SQL in a file
    local $/ = undef;
    open(SQLTMP, $f_sql) || die;
    $sql = <SQLTMP>;
    close SQLTMP;
  } else {
    confess "specify -sql or -table";
  }
  die "specify -table or -outfile" unless $outfile;
  my $query = $dbi->prepare($sql) || die;
  $query->execute() || die;

  # don't use Reporter.pm, because it's possible the query might 
  # return multiple columns with the same label

  my $delimiter = "\t";
  open(EQ, ">" . $outfile) || die;
  printf EQ "%s\n", join $delimiter, @{$query->{NAME}} if $header;
  while (my $row = $query->fetchrow_arrayref) {
    my @row = @{$row};
    if ($trim_whitespace) {
      # strip any leading and trailing whitespace
      foreach (@row) {
	s/^\s+//;
	s/\s+$//;
      }
    }
    printf EQ "%s\n", join $delimiter, @row;
  }
  close EQ;
}

sub parse_pg_pass {
  # parse PostreSQL ~/.pgpass file
  my $f_pgpass = $ENV{HOME} . "/.pgpass";
  die "where is $f_pgpass" unless -s $f_pgpass;
  # TO DO: cleanup exceptions

  open(PGPTMP, $f_pgpass) || die;
  my $hdr = <PGPTMP>;
  chomp $hdr;
  die unless $hdr =~ /^#/;
  $hdr =~ s/^#//;
  my @h = split /:/, $hdr;
  die unless @h == 5;
  # ever any variation?

  my %host2info;;
  while (my $line = <PGPTMP>) {
    chomp $line;
    my @f = split /:/, $line;
    die unless @f == @h;
    my %r;
    @r{@h} = @f;

    my $host = $r{hostname} || die;
    die "duplicate" if $host2info{$host};
    $host2info{$host} = \%r;

  }
  close PGPTMP;
  return \%host2info;
}

sub parse_pg_service {
  # parse PostreSQL ~/.pg_service.conf file
  # - also patching in passwords from pgpass
  my %options = @_;
  # STATIC
  my $f_pg_service = $ENV{HOME} . "/.pg_service.conf";
  die "where is $f_pg_service" unless -s $f_pg_service;
  # TO DO: cleanup exceptions

  my %services;

  open(PGSTMP, $f_pg_service) || die;
  my $service_name;
  while (<PGSTMP>) {
    chomp;

    if (/^\[(.*)\]$/) {
      $service_name = $1;
    } elsif (/\w/) {
      die unless $service_name;
      my @f = split /=/, $_;
      die unless @f == 2;
      $services{$service_name}{$f[0]} = $f[1];
      # key/value pair
    }
  }

  my $pgpass = parse_pg_pass();

  foreach my $sname (sort keys %services) {
    my $service = $services{$sname} || die;
    my $s_host = $service->{host} || die;
    my $s_port = $service->{port} || die;
    my $s_dbname = $service->{dbname} || die;
    my $s_user = $service->{user} || die;

    if (my $pp = $pgpass->{$s_host}) {
      # .pgpass entry found for this host
      if ($pp->{database} eq $s_dbname and
	  $pp->{port} == $s_port and
	  $pp->{username} eq $s_user) {
	# .pgpass details match service details, so usable
	$service->{password} = $pp->{password};
#	printf STDERR "set service password for %s\n", $sname;
      }
    }
  }

  my $result;
  if (my $sname = $options{"-get-dbi"}) {
    # create DBI connection to specified service
    my $sref = $services{$sname} || die "no service $sname";
#    dump_die($sref, "debug", 1);
    my @dsn;
    foreach my $k (qw(dbname host port)) {
      push @dsn, sprintf '%s=%s', $k, $sref->{$k} || die;
    }
    my $dsn = sprintf "dbi:Pg:%s", join ";", @dsn;
    my $user = $sref->{user} || die;
    my $pw = $sref->{password} || die;
    $result = DBI->connect($dsn, $user, $pw) || die "db connect failed";
  } else {
    $result = \%services;
  }
  return $result;
}

1;
