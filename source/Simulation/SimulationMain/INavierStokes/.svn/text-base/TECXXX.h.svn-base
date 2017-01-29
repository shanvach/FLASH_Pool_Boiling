/*
 * TECXXX.h: Copyright (C) 1988-2003 Amtec Engineering, Inc.
 */

#if !defined TECXXX_H_
#define TECXXX_H_

#if !defined CRAY
#  define TECFOREIGN100 tecforeign100
#  define TECINI100     tecini100
#  define TECZNE100     teczne100
#  define TECDAT100     tecdat100
#  define TECNOD100     tecnod100
#  define TECGEO100     tecgeo100
#  define TECTXT100     tectxt100
#  define TECLAB100     teclab100
#  define TECFIL100     tecfil100
#  define TECEND100     tecend100
#  define TECUSR100     tecusr100
#  define TECAUXSTR100  tecauxstr100
#  define TECZAUXSTR100 teczauxstr100
#  define TECVAUXSTR100 tecvauxstr100
#  define TECFACE100    tecface100

#  define TECINI  tecini
#  define TECZNE  teczne
#  define TECDAT  tecdat
#  define TECNOD  tecnod
#  define TECGEO  tecgeo
#  define TECTXT  tectxt
#  define TECLAB  teclab
#  define TECFIL  tecfil
#  define TECEND  tecend
#  define TECUSR  tecusr
#endif



#define INTEGER4  int

#if !defined MSWIN
# if defined (_WINDOWS) || defined (WIN32)
#   define MSWIN
# endif /* _WINDOWS || WIN32 */
#endif /* !MSWIN */

#if !defined (EXTERNC)
# if defined (__cplusplus)
#  define EXTERNC extern "C"
# else
#  define EXTERNC
# endif /* __cplusplus */
#endif /* EXTERN_C */

#if !defined (STDCALL)
# if defined MSWIN
#  define STDCALL __stdcall
# else /* !MSWIN */
#  define STDCALL
# endif /* MSWIN */
#endif /* STDCALL */

#if !defined (DLLEXPORT)
# if defined (MSWIN)
#  define DLLEXPORT _declspec (dllexport)
# else
#  define DLLEXPORT
# endif /* MSWIN */
#endif /* DLLEXPORT */

#if !defined (DLLIMPORT)
# if defined (MSWIN)
#  define DLLIMPORT _declspec (dllimport)
# else
#  define DLLIMPORT
# endif /* MSWIN */
#endif /* DLLIMPORT */


#if defined (TECPLOTKERNEL)
/* CORE SOURCE CODE REMOVED */
#else /* !TECPLOTKERNAL && !MAKEARCHIVE */
# define LIBCALL STDCALL
# define LIBFUNCTION EXTERNC DLLIMPORT
#endif

/* Old V9 functions retained for backward compatibility */

LIBFUNCTION INTEGER4 LIBCALL TECINI(char     *Title,
                                    char     *Variables,
                                    char     *FName,
                                    char     *ScratchDir,
                                    INTEGER4 *Debug,
                                    INTEGER4 *VIsDouble);

LIBFUNCTION INTEGER4 LIBCALL TECZNE(char     *ZoneTitle,
                                    INTEGER4 *IMx,
                                    INTEGER4 *JMx,
                                    INTEGER4 *KMx,
                                    char     *ZFormat,
                                    char     *DupList);

LIBFUNCTION INTEGER4 LIBCALL TECDAT(INTEGER4  *N,
                                    void      *FieldData,
                                    INTEGER4  *IsDouble);

LIBFUNCTION INTEGER4 LIBCALL TECNOD(INTEGER4 *NData);

LIBFUNCTION INTEGER4 LIBCALL TECEND(void);

LIBFUNCTION INTEGER4 LIBCALL TECLAB(char *S);

LIBFUNCTION INTEGER4 LIBCALL TECUSR(char *S);

LIBFUNCTION INTEGER4 LIBCALL TECGEO(double    *XPos,
                                    double    *YPos,
                                    double    *ZPos,
                                    INTEGER4  *PosCoordMode,
                                    INTEGER4  *AttachToZone,
                                    INTEGER4  *Zone,
                                    INTEGER4  *Color,
                                    INTEGER4  *FillColor,
                                    INTEGER4  *IsFilled,
                                    INTEGER4  *GeomType,
                                    INTEGER4  *LinePattern,
                                    double    *PatternLength,
                                    double    *LineThickness,
                                    INTEGER4  *NumEllipsePts,
                                    INTEGER4  *ArrowheadStyle,
                                    INTEGER4  *ArrowheadAttachment,
                                    double    *ArrowheadSize,
                                    double    *ArrowheadAngle,
                                    INTEGER4  *Scope,
                                    INTEGER4  *NumSegments,
                                    INTEGER4  *NumSegPts,
                                    float     *XGeomData,
                                    float     *YGeomData,
                                    float     *ZGeomData,
                                    char      *mfc);

LIBFUNCTION INTEGER4 LIBCALL TECTXT(double    *XPos,
                                    double    *YPos,
                                    INTEGER4  *PosCoordMode,
                                    INTEGER4  *AttachToZone,
                                    INTEGER4  *Zone,
                                    INTEGER4  *BFont,
                                    INTEGER4  *FontHeightUnits,
                                    double    *FontHeight,
                                    INTEGER4  *BoxType,
                                    double    *BoxMargin,
                                    double    *BoxLineThickness,
                                    INTEGER4  *BoxColor,
                                    INTEGER4  *BoxFillColor,
                                    double    *Angle,
                                    INTEGER4  *Anchor,
                                    double    *LineSpacing,
                                    INTEGER4  *TextColor,
                                    INTEGER4  *Scope,
                                    char      *Text,
                                    char      *mfc);

LIBFUNCTION INTEGER4 LIBCALL TECFIL(INTEGER4 *F);


/* New V10 tecio functions */


LIBFUNCTION void LIBCALL TECFOREIGN100(INTEGER4 *OutputForeignByteOrder);

LIBFUNCTION INTEGER4 LIBCALL TECINI100(char     *Title,
                                       char     *Variables,
                                       char     *FName,
                                       char     *ScratchDir,
                                       INTEGER4 *Debug,
                                       INTEGER4 *VIsDouble);

LIBFUNCTION INTEGER4 LIBCALL TECZNE100(char     *ZoneTitle,
                                       INTEGER4 *ZoneType,
                                       INTEGER4 *IMxOrNumPts,
                                       INTEGER4 *JMxOrNumElements,
                                       INTEGER4 *KMx,
                                       INTEGER4 *ICellMx,
                                       INTEGER4 *JCellMx,
                                       INTEGER4 *KCellMx,
                                       INTEGER4 *IsBlock,
                                       INTEGER4 *NumFaceConnections,
                                       INTEGER4 *FaceNeighborMode,
                                       INTEGER4 *ValueLocation,
                                       INTEGER4 *ShareVarFromZone,
                                       INTEGER4 *ShareConnectivityFromZone);

LIBFUNCTION INTEGER4 LIBCALL TECDAT100(INTEGER4  *N,
                                       void      *FieldData,
                                       INTEGER4  *IsDouble);

LIBFUNCTION INTEGER4 LIBCALL TECNOD100(INTEGER4 *NData);

LIBFUNCTION INTEGER4 LIBCALL TECEND100(void);

LIBFUNCTION INTEGER4 LIBCALL TECLAB100(char *S);

LIBFUNCTION INTEGER4 LIBCALL TECUSR100(char *S);

LIBFUNCTION INTEGER4 LIBCALL TECGEO100(double    *XPos,
                                       double    *YPos,
                                       double    *ZPos,
                                       INTEGER4  *PosCoordMode,
                                       INTEGER4  *AttachToZone,
                                       INTEGER4  *Zone,
                                       INTEGER4  *Color,
                                       INTEGER4  *FillColor,
                                       INTEGER4  *IsFilled,
                                       INTEGER4  *GeomType,
                                       INTEGER4  *LinePattern,
                                       double    *PatternLength,
                                       double    *LineThickness,
                                       INTEGER4  *NumEllipsePts,
                                       INTEGER4  *ArrowheadStyle,
                                       INTEGER4  *ArrowheadAttachment,
                                       double    *ArrowheadSize,
                                       double    *ArrowheadAngle,
                                       INTEGER4  *Scope,
                                       INTEGER4  *Clipping,
                                       INTEGER4  *NumSegments,
                                       INTEGER4  *NumSegPts,
                                       float     *XGeomData,
                                       float     *YGeomData,
                                       float     *ZGeomData,
                                       char      *mfc);

LIBFUNCTION INTEGER4 LIBCALL TECTXT100(double    *XOrThetaPos,
                                       double    *YOrRPos,
                                       double    *ZOrUnusedPos,
                                       INTEGER4  *PosCoordMode,
                                       INTEGER4  *AttachToZone,
                                       INTEGER4  *Zone,
                                       INTEGER4  *BFont,
                                       INTEGER4  *FontHeightUnits,
                                       double    *FontHeight,
                                       INTEGER4  *BoxType,
                                       double    *BoxMargin,
                                       double    *BoxLineThickness,
                                       INTEGER4  *BoxColor,
                                       INTEGER4  *BoxFillColor,
                                       double    *Angle,
                                       INTEGER4  *Anchor,
                                       double    *LineSpacing,
                                       INTEGER4  *TextColor,
                                       INTEGER4  *Scope,
                                       INTEGER4  *Clipping,
                                       char      *String,
                                       char      *mfc);

LIBFUNCTION INTEGER4 LIBCALL TECFIL100(INTEGER4 *F);

LIBFUNCTION INTEGER4 LIBCALL TECAUXSTR100(char *Name,
                                          char *Value);

LIBFUNCTION INTEGER4 LIBCALL TECZAUXSTR100(char *Name,
                                           char *Value);

LIBFUNCTION INTEGER4 LIBCALL TECVAUXSTR100(INTEGER4 *Var,
                                           char     *Name,
                                           char     *Value);

LIBFUNCTION INTEGER4 LIBCALL TECFACE100(INTEGER4 *FaceConnections);

#endif /* TECXXX_H_ */
