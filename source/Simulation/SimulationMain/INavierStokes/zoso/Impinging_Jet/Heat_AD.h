#define DIR_X 1
#define DIR_Y 2
#define DIR_Z 3

#define SCHEME 1

#if SCHEME == 1
#define HEAT_CENTRAL
#endif

#if SCHEME == 2
#define HEAT_UPWIND
#endif
