#! /usr/bin/perl
BEGIN{$lastEle = ''; $iprop = 0;}
LINE:
    while (<>) {
	chomp;
#	print "NOW HANDLING \"$_\"...\n";
	$propertyNameIn = $_;		# particle property name
	$unkName = $_;		# UNK name
	if (m/^(..).+(..)$/) {
	    $unkName = "$1$2";
	}
	if (m/^([A-Za-z]{1,2})([12][0-9][0-9]|[1-9][0-9]*)/) {
	    $isAbund = 1;
	    $unkEle = $1;
	} elsif (m/^n$/) {
	    $isAbund = 1;
	    $unkEle = "neu";
	} else {
	    $isAbund = 0;
	    $unkEle = '';
	}


	$preDefinedInConfigTop = 0;
	$typeInfo = ""; 		# type info
	if (m/(dens|pres|temp|eint|vel[xyz]|entr)/) {
	    $variableKind = "VARIABLE"; # variable kind
	    $isAbund = 0;
	    if (m/dens/) { $typeInfo = " TYPE: PER_VOLUME"; }
	    elsif (m/(eint|vel[xyz]|entr)/) { $typeInfo = " TYPE: PER_MASS"; }
	    if (m/(dens|temp|eint)/) {
		$preDefinedInConfigTop = 1;
	    }
	}
	elsif (m/(sumy|ye)/) {
	    $variableKind = "MASS_SCALAR"; # variable kind
	    $isAbund = 0;
	    $preDefinedInConfigTop = 1;
	}
	elsif (m/(rpv[0-9])/) {$variableKind = "MASS_SCALAR";}
	elsif (m/(flam)/) {$variableKind = "MASS_SCALAR";}
	else {$variableKind = "SPECIES";}
#	$iprop = $.; 

#	print "UNKELE=$unkEle, LASTELE=$lastEle\n";
	if ($isAbund) {
#	    if ($unkEle ne $lastEle) {
	    if (! defined $speciesDefined{$unkEle}) {
		print $variableKind . " " . $unkEle . $typeInfo . "\n";
	    }
	} else {
	    if (! $preDefinedInConfigTop) {
		print $variableKind . " " . $unkName . $typeInfo . "\n";
	    } else {
		print "####" . $variableKind . " " . $unkName . $typeInfo . "\n";
		if ($variableKind eq "MASS_SCALAR") {
		    print "\n";
		    next LINE;
		    }
	    }
	}
	$iprop++; 
	print "PARTICLEPROP " . $propertyNameIn . " REAL\n";
	if ($isAbund) {
	    $unkNameOut = $unkEle;
	} else {
	    $unkNameOut = $unkName;
	}
	print "PARTICLEMAP TO " . $propertyNameIn . " FROM " . $variableKind . " $unkNameOut\n";
	print "PARAMETER particle_attribute_$iprop STRING \"\L$propertyNameIn\"\n";
	print "\n";
	if ($isAbund) {$lastEle = $unkEle;}
	if (! defined $speciesDefined{$unkEle}) {$speciesDefined{$unkEle} = 1;}
    }
