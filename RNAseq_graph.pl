#!/usr/bin/perl -w
use strict;

## genes.fpkm_tracking
## tracking_id	class_code	nearest_ref_id	gene_id	gene_short_name	tss_id	locus	length	coverage	q1_FPKM	q1_conf_lo	q1_conf_hi	q1_status	q2_FPKM	q2_conf_lo	q2_conf_hi	q2_status	q3_FPKM	q3_conf_lo	q3_conf_hi	q3_status
## ChrSy.fgenesh.gene.1	-	-	ChrSy.fgenesh.gene.1	ChrSy.fgenesh.gene.1	-	chrSy:784-1168	-	-	0	0	0	OK	0	0	0	OK	0	0	0	OK

my $top_dir = shift;
opendir my ($dh), $top_dir or die "Couldn't open dir '$top_dir': $!";
my @files = readdir $dh;
closedir $dh;

my %order;
my %fpkm;
foreach my $dir (@files) {
  next unless -d "$top_dir/$dir";
  my ( $run_info, $fpkm_file );
  if ( -e "$top_dir/$dir/run.info" and -e "$top_dir/$dir/genes.fpkm_tracking" )
  {
    $run_info  = "$top_dir/$dir/run.info";
    $fpkm_file = "$top_dir/$dir/genes.fpkm_tracking";
  } else {
    next;
  }
  my @order;
  my $group;
  open INFO, $run_info or die "Can't open $run_info\n";
  while ( my $line = <INFO> ) {
    chomp $line;
    next unless $line =~ /cmd_line/;
    ($group) = $line =~ /\-o\s+(\S+)/;
    my @bams = $line =~ /(\S+\.bam)/g;
    foreach my $bam (@bams) {
      if ( $bam =~ /NB/ ) {
        push @order, 'NB';
      } elsif ( $bam =~ /HEG/ ) {
        push @order, 'HEG4';
      } else {
        push @order, 'EG4';
      }
    }
  }
  $order{$group}{ $order[0] } = 'q1';
  $order{$group}{ $order[1] } = 'q2';
  $order{$group}{ $order[2] } = 'q3';

  open IN, $fpkm_file or die "Can't open $fpkm_file\n";
  while ( my $line = <IN> ) {
    chomp $line;
    next if $line =~ /^tracking_id/;
    my @line    = split "\t", $line;
    my $gene    = $line[0];
    my $q1_FPKM = $line[9];
    my $q2_FPKM = $line[13];
    my $q3_FPKM = $line[17];
    foreach my $strain ( keys %{ $order{$group} } ) {

      if ( $order{$group}{$strain} eq 'q1' ) {
        $fpkm{$gene}{$group}{$strain} = $q1_FPKM;
      } elsif ( $order{$group}{$strain} eq 'q2' ) {
        $fpkm{$gene}{$group}{$strain} = $q2_FPKM;
      } elsif ( $order{$group}{$strain} eq 'q3' ) {
        $fpkm{$gene}{$group}{$strain} = $q3_FPKM;
      }
    }
  }
}
my $count     = 0;
#my $dir_count = 0;
my $new_dir   = "RNASeq_images";
if ( !-d $new_dir ) {
  `mkdir -p $new_dir`;
}
foreach my $gene ( sort keys %fpkm ) {
  my @dataframe;
  my ($subdir) = $gene =~ /^(.{8})/; #LOC_Os1 
  my $new_dir = "RNASeq_images/$subdir";
  my @values;
  `mkdir -p $new_dir` if !-e $new_dir;
  next if -e "$new_dir/$gene.RNASeq.png";
  open OUTR, ">$new_dir/$gene.RNASeq.R"
    or die "Cant open $new_dir/$gene.RNASeq.R";

  print OUTR 
qq(png(file = "$new_dir/$gene.RNASeq.png", bg = "transparent")\n
par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)\n);

  #my %dataframe;
  foreach my $group ( sort keys %{ $fpkm{$gene} } ) {
    my $NB_FPKM   = $fpkm{$gene}{$group}{NB};
    my $HEG4_FPKM = $fpkm{$gene}{$group}{HEG4};
    my $EG4_FPKM  = $fpkm{$gene}{$group}{EG4};
    
    print OUTR "$group <- c(", join( ',', $NB_FPKM, $EG4_FPKM, $HEG4_FPKM ),")\n";

    push @dataframe, $group;
  }
  my $colors = scalar @dataframe;
  my $dataframe = join ",", @dataframe;

  print OUTR
qq(colorlist <- append ( colors()[grep("medium",colors())]  , colors()[grep("dark",colors())] )\n
plot_colors <- c(colorlist[1:$colors])\n
data <- data.frame($dataframe)\n
max_y <- max(data)\n);

  my $first_dataframe = shift @dataframe;
  my $pch             = 0;

  print OUTR
qq(plot(data\$$first_dataframe,type = "b",pch=$pch,col=plot_colors[1],ylim=c(0,max_y),axes=FALSE,ann=FALSE)\n
axis(1,at=1:3,lab=c(\"NB\",\"EG4\",\"HEG4\"))\n
axis(2,las=1, at=1*0:max_y)\n
box()\n);

  my $i = 2;

  foreach my $df (@dataframe) {
    if ( $pch == 25 ) {
      $pch = 0;
    } else {
      $pch++;
    }
    print OUTR qq(points(data\$$df, type="b",pch=$pch,col=plot_colors[$i])\n);
    $i++;
  }
  print OUTR qq ( 
title(main= paste(strwrap("Gene expression in different Conditions and Strains for $gene",width=50),collapse="\\n"),col.main="black", font.main=2)\n
title(xlab= "Strains", col.lab=rgb(0,0,0))\n
title(ylab= "FPKM", col.lab=rgb(0,0,0))\n
legend(3.1, max_y, names(data), cex=0.8, col=plot_colors, pch=0:$colors)\n
dev.off()\n);
  close OUTR;

  `Rscript "$new_dir/$gene.RNASeq.R"`;
}

__END__
#my %conditions =
#(cold_day1_1hr_control => "aliceblue", "cold_day1_3hr_control" => "azure","cold_day2_1hr_control"=>"blue","cold_day2_3hr_control"=>"blueviolet","cold_day3_1hr_control"=>"cadetblue","cold_day3_3hr_control"=>"cornflowerblue","cold_day1_1hr_exp"=>"cyan","cold_day1_3hr_exp"=>"darkblue","cold_day2_1hr_exp"=>"darkcyan","cold_day2_3hr_exp"=>"darkslate","cold_day3_1hr_exp"=>"darkturquois","cold_day3_3hr_exp"=>"deepskyblue","salt_0.5hr_control"=>"antiquewhite","salt_2hr_control"=>"beige","salt_10hr_control"=>"bisque","salt_24hr_control"=>"blanchedalmond","salt_48hr_control"=>"burlywood","salt_30min_control"=>"burlywood1","salt_30min_150mM"=>"burlywood4","salt_2hr_150mM"=>"wheat","salt_10hr_150mM"=>"floralwhite","salt_24hr_150mM"=>"lemonchiffon","salt_48hr_150mM"=>"lightgoldenrod","drought_0hr_control"=>"violetred1","drought_1hr_control"=>"tomato","drought_5hr_exp"=>"pink","drought_10hr_exp"=>"paleviletred","drought_24hr_exp"=>"salmon","drought_48hr_exp"=>"rosybrown1","drought_a"=>"orangered","drought_b"=>"red","drought_c"=>"lightsalmon","drought_d"=>"lightpink","drought_e"=>"lightcoral","drought_f"=>"indianred","drought_g"=>"hotpink");
