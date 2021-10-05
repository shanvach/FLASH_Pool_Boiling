#!MC 1410
$!VarSet |MFBD| = '/nobackup/adhruv/FLOW/IOData'

$!VarSet |FileIndex| = 4

$! LOOP 1

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

$! LOOP 1399

      $!VarSet |proc| = |LOOP|

      $!IF |proc| < 10
      $!VarSet |procIndex| = "000|proc|"
      $!ELSEIF |proc| < 100
      $!VarSet |procIndex| = "00|proc|"
      $!ELSEIF |proc| < 1000
      $!VarSet |procIndex| = "0|proc|"
      $!ELSEIF |proc| < 10000
      $!VarSet |procIndex| = "|proc|"
      $!ENDIF

$!READDATASET  '"|MFBD|/data.|fINDEX|.|procIndex|.plt"'
  READDATAOPTION = APPEND
  RESETSTYLE = YES
  VARLOADMODE = BYNAME
  ASSIGNSTRANDIDS = YES

$!ENDLOOP

$!VarSet |fIndex| +=14

$!READSTYLESHEET  "|MFBD|/flow_boiling_iso_temp.sty"
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
$!EXPORTSETUP EXPORTFNAME = '|MFBD|/images/ISO-TEMP|fIndex|.png'
$!EXPORT 
  EXPORTREGION = CURRENTFRAME

$!READSTYLESHEET  "|MFBD|/flow_boiling_iso_bub.sty"
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
$!EXPORTSETUP EXPORTFNAME = '|MFBD|/images/ISO-BUB|fIndex|.png'
$!EXPORT 
  EXPORTREGION = CURRENTFRAME

$!ENDLOOP

$!RemoveVar |MFBD|
