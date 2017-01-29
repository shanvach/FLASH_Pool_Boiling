#!/usr/bin/perl -w

# Convergence study for the Noh problem

my $title    = "Noh";
my $term     = "png";
my $dotterm  = ".png";
my $lw       = 2;

my $dumps_dir = "dumps_ppm";

my $dekfile  = "flash.par";
my $parfile  = "this.par";

my $idl      = "/usr/local/rsi/idl_6.1/bin/idl";   # Path to idl

my $r_max = 1.0;

my $max_dumps = 60;                                # Last checkpoint file number
my @colors = (5,4,2,3,1);                          # gnuplot line styles
my @n_blocks = (100,250,500,750,1000);             # Number of blocks


my %sim_params;

my @r=(); 
my @rho=();
my @p=(); 
my @u=();

my @rho_sol=();
my @p_sol=();
my @u_sol=();

my $nb;


&clean_start;
foreach $nb (@n_blocks) {
    &clean_dir;
    &set_sim_params($nb);
    &make_dek($dekfile,$parfile);
     
    $status = system "nohup", "mpirun", "-np", "2", "flash3", "-par_file", $parfile, "&";
    die if ($status);
  
    $status = system $idl, "idlrunlast";
    die if ($status);

    &process_data($nb,$max_dumps);
}
&make_plots;

die if system "mv", "dumps", $dumps_dir;
@cplist = glob "saved/*";
foreach (@cplist) {
    die if system "cp", $_, "$dumps_dir/";
}



sub clean_start {
    die if system "rm", "-rf", $dumps_dir;
    die if system "rm", "-rf", "dumps";
    die if system "mkdir", "dumps";
}


sub clean_dir {
    my @dlist = glob "dumps/noh_* dumps/rho_* dumps/p_* dumps/u_* dumps/rho-* dumps/p-* dumps/u-*";
    foreach (@dlist) {
	unlink $_;
    }
}


sub set_sim_params {
  my $nb = shift;

  $sim_param{"nblocks"} = $nb;
}


sub make_dek {
  my $sfilename = shift;
  my $dfilename = shift;

  open INDEK, "<$sfilename" or die;
  open OUTDEK, ">$dfilename" or die;

  while (<INDEK>) {
      chomp;

      if (/nblockx/i) {
	  print OUTDEK "nblockx = ", $sim_param{"nblocks"} . "\n";
      }
#      elsif (/nblocky/i) {
#	  print OUTDEK "nblocky = ", $sim_param{"nblocks"} . "\n";	  
#      }
#      elsif (/nblockz/i) {
#	  print OUTDEK "nblockz = ", $sim_param{"nblocks"} . "\n";
#      }
      else {
	  print OUTDEK "$_\n";
      }
  }

  close INDEK;
  close OUTDEK;
}





sub collect_dumps {
    my $dnum = shift;
    my $n;
    my $t;
    my $i;
    my $x;
    my $a;

    @r=(); 
    @rho=();
    @p=(); 
    @u=();
    
    $i = 0;

    print "dnum = $dnum\n";
    open RHOFILE, "<dumps/rho_$dnum" or die;
    while (<RHOFILE>) {
	chomp;
	if (!/\#/i) {
	    ($x,$a)  = split(" ",$_);
            if ($x <= $r_max) {
		$r[$i]   = $x;
		$rho[$i] = $a;
		$i++; 
	    }
	}
	else {
	    ($s,$n,$t) = split(" ",$_);
	}
    }
    close RHOFILE;

    $i = 0;
    open PFILE, "<dumps/p_$dnum" or die;
    while (<PFILE>) {
	chomp;
	if (!/\#/i) {
	    ($x,$a)  = split(" ",$_);
	    if ($x <= $r_max) {
		$p[$i]   = $a;
		$i++; 
	    }
	}
    }
    close PFILE;

    $i = 0;
    open UFILE, "<dumps/u_$dnum" or die;
    while (<UFILE>) {
	chomp;
	if (!/\#/i) {
	    ($x,$a)  = split(" ",$_);
	    if ($x <= $r_max) {
		$u[$i]   = $a;
		$i++; 
	    }
	}
    }
    close UFILE;

    $n = @r;

    open RFILE, ">dumps/r.dat" or die;
    print RFILE "# $n\n";
    for ($i=0; $i<$n; $i++) {
	print RFILE "$r[$i]\n";
    }
    close RFILE;

    ($n,$t);
}



sub make_sol {
    my $t = shift;
    my $filename = shift;
 
    @rho_sol = ();
    @p_sol = ();
    @u_sol = ();
  
    my $i = 0;
    my $n = @r;

    open SOLFILE, ">$filename" or die;
    
    for ($i=0; $i<$n; $i++) {
	my $dens;
	my $pres;
	my $velx;

	if ($r[$i] <= 0.2) {
	    $dens = 64;
	    $pres = 21.3333333333;
            $velx = 0;
	}
	else {
	    $dens = (1.0 + 0.6/$r[$i])**2;
	    $pres = 0;
	    $velx = -1;
	}

	$rho_sol[$i] = $dens;
	$p_sol[$i]   = $pres;
	$u_sol[$i]   = $velx;
            
	print SOLFILE "$r[$i] $rho_sol[$i] $p_sol[$i] $u_sol[$i]\n";
    }
    close SOLFILE;
}



sub calc_errors {
    my $nb = shift;
    my $filename = shift;

    my $rho_sum = 0.0;
    my $p_sum = 0.0;
    my $u_sum = 0.0;

    my $srho = 0.0;
    my $sp = 0.0;
    my $su = 0.0;
    
    my $i;
    my $n;

    my $r1;
    my $r2;
    my $dr;

    $n = @r;
    for ($i=0; $i<($n-1); $i++) {
	if ($i==0) {
	    $r1 = 0.0;
	}
	$r2 = 0.5*($r[$i+1] + $r[$i]);
	$dr = $r2 - $r1;

	$rho_sum = $rho_sum + abs($rho[$i] - $rho_sol[$i])*$dr;
	$p_sum = $p_sum + abs($p[$i] - $p_sol[$i])*$dr;
	$u_sum = $u_sum + abs($u[$i] - $u_sol[$i])*$dr;

	$srho = $srho + $rho_sol[$i]*$dr;
	$sp = $sp +  $p_sol[$i]*$dr;
	$su = $su + $u_sol[$i]*$dr;

	$r1 = $r2;
    }

#    print "u_sum = $u_sum \n";
#    print "su = $su \n";
    
    $rho_sum = $rho_sum / $srho;
    $p_sum = $p_sum / $sp;
    $u_sum = $u_sum / abs($su);

#    print "norm = $u_sum \n";

    open EFILE, ">>$filename" or die;
    print EFILE "$nb $rho_sum $p_sum $u_sum\n";
    close EFILE;
}



sub process_data {
    my $nb = shift;
    my $max_dnum = shift;

    my $n;
    my $t;

    my $solfile = "dumps/nohsol.dat";
    my $errfile = "dumps/error.dat";

    ($n,$t) = &collect_dumps($max_dnum);
	
    &make_sol($t,$solfile);
    &calc_errors($nb,$errfile);

    system "cp", "dumps/rho_$max_dnum", "dumps/rhonb" . $nb;
    system "cp", "dumps/p_$max_dnum", "dumps/pnb" . $nb;
    system "cp", "dumps/u_$max_dnum", "dumps/unb" . $nb;    
}



sub make_plots {
    my $solfile  = "dumps/nohsol.dat";

    open PLOTFILE, ">dumps/plot_data" or die;
    print PLOTFILE "set term $term\n";
    print PLOTFILE "set key top right\n";
    print PLOTFILE "set key box\n";
    print PLOTFILE "set output \"dumps/$title" . "_rho" . "$dotterm\"\n";
    print PLOTFILE "set title \"$title: Density\"\n";
    print PLOTFILE "set xlabel \"r\"\n";
    print PLOTFILE "set ylabel \"rho\"\n";
    print PLOTFILE "set xrange [0:1]\n";
    print PLOTFILE "plot \"$solfile\" using 1:2 w l lw $lw lt -1 t \"Solution\"";
    for ($i=0; $i<@n_blocks; $i++) {
	$dfile = "dumps/rhonb" . $n_blocks[$i];
	print PLOTFILE ", \"$dfile\" using 1:2 w l lw $lw lt $colors[$i]  t \"$n_blocks[$i] blocks\"";
    }
    close PLOTFILE;    
    system "gnuplot", "dumps/plot_data";
    unlink "dumps/plot_data";

    open PLOTFILE, ">dumps/plot_data" or die;
    print PLOTFILE "set term $term\n";
    print PLOTFILE "set key top right\n";
    print PLOTFILE "set key box\n";
    print PLOTFILE "set output \"dumps/$title" . "_p" . "$dotterm\"\n";
    print PLOTFILE "set title \"$title: Pressure\"\n";
    print PLOTFILE "set xlabel \"r\"\n";
    print PLOTFILE "set ylabel \"p\"\n";
    print PLOTFILE "set xrange [0:1]\n";
    print PLOTFILE "plot \"$solfile\" using 1:3 w l lw $lw lt -1 t \"Solution\"";
    for ($i=0; $i<@n_blocks; $i++) {
	$dfile = "dumps/pnb" . $n_blocks[$i];
	print PLOTFILE ", \"$dfile\" using 1:2 w l lw $lw lt $colors[$i]  t \"$n_blocks[$i] blocks\"";
    }
    close PLOTFILE;    
    system "gnuplot", "dumps/plot_data";
    unlink "dumps/plot_data";

   open PLOTFILE, ">dumps/plot_data" or die;
    print PLOTFILE "set term $term\n";
    print PLOTFILE "set key top right\n";
    print PLOTFILE "set key box\n";
    print PLOTFILE "set output \"dumps/$title" . "_u" . "$dotterm\"\n";
    print PLOTFILE "set title \"$title: Velocity\"\n";
    print PLOTFILE "set xlabel \"r\"\n";
    print PLOTFILE "set ylabel \"u\"\n";
    print PLOTFILE "set xrange [0:1]\n";
    print PLOTFILE "plot \"$solfile\" using 1:4 w l lw $lw lt -1 t \"Solution\"";
    for ($i=0; $i<@n_blocks; $i++) {
	$dfile = "dumps/unb" . $n_blocks[$i];
	print PLOTFILE ", \"$dfile\" using 1:2 w l lw $lw lt $colors[$i]  t \"$n_blocks[$i] blocks\"";
    }
    close PLOTFILE;    
    system "gnuplot", "dumps/plot_data";
    unlink "dumps/plot_data";
}
