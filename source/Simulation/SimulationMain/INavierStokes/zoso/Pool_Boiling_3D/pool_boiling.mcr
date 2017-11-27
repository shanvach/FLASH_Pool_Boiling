#!MC 1410
$!VarSet |MFBD| = '/lustre/groups/balarasgrp/akash/sim_NB_Earth/IOData_Fr1_Tw80_new_version1'

$!VarSet |FileIndex| = -1

$! LOOP 366

        $!VarSet |FileIndex| += 1

        $!IF |FileIndex| < 10
        $!VarSet |fIndex| = "000|FileIndex|"
        $!ELSEIF |FileIndex| < 100
        $!VarSet |fIndex| = "00|FileIndex|"
        $!ELSEIF |FileIndex| < 1000
        $!VarSet |fIndex| = "0|FileIndex|"
        $!ELSEIF  |FileIndex| < 10000
        $!VarSet |fIndex| = "|FileIndex|"
        $!ENDIF

$!READDATASET  '"|MFBD|/data.|fINDEX|.0000.plt"'
  READDATAOPTION = NEW
  RESETSTYLE = YES
  VARLOADMODE = BYNAME
  ASSIGNSTRANDIDS = YES

$! LOOP 143

      $!VarSet |proc| = |LOOP|

      $!IF |proc| < 10
      $!VarSet |procIndex| = "000|proc|"
      $!ELSEIF |proc| < 100
      $!VarSet |procIndex| = "00|proc|"
      $!ELSEIF |proc| < 200
      $!VarSet |procIndex| = "0|proc|"
      $!ENDIF

$!READDATASET  '"|MFBD|/data.|fINDEX|.|procIndex|.plt"'
  READDATAOPTION = APPEND
  RESETSTYLE = YES
  VARLOADMODE = BYNAME
  ASSIGNSTRANDIDS = YES

$!ENDLOOP

$!VarSet |fIndex| +=0

$!READSTYLESHEET  "|MFBD|/pool_boiling.sty"
  INCLUDEPLOTSTYLE = YES
  INCLUDETEXT = YES
  INCLUDEGEOM = YES
  INCLUDEAUXDATA = YES
  INCLUDESTREAMPOSITIONS = YES
  INCLUDECONTOURLEVELS = YES
  MERGE = NO
  INCLUDEFRAMESIZEANDPOSITION = NO
$!PRINTSETUP PALETTE = COLOR
$!EXPORTSETUP IMAGEWIDTH = 3000
$!EXPORTSETUP EXPORTFNAME = '|MFBD|/images/BUB|fIndex|.png'
$!EXPORT 
  EXPORTREGION = CURRENTFRAME

$!ENDLOOP

$!RemoveVar |MFBD|
