
/*----- CPP defines for coral module ----------------*/
/***  Biological model options. ***/

#define REEF_ECOSYS  /*Original CPP flag */

#if defined REEF_ECOSYS
# define BIOLOGY
# define ANA_BIOLOGY

/* compartments */
# define ORGANIC_MATTER  /*Original CPP flag */
/*# define CARBON_ISOTOPE*/  /*Original CPP flag */
# define NUTRIENTS   /*Original CPP flag */

/*# define CORAL_POLYP*/  /* USE coral module */
/*# define SEAGRASS*/     /* USE seagrass module */
/*# define ALGAE*/        /* USE algae module  */
# define SEDIMENT_ECOSYS        /* USE sedecosys module  */

# if defined ORGANIC_MATTER
/*#  define FOODWEB*/      /* USE foodweb module */
# endif
# if defined FOODWEB
#  define COT_STARFISH      /* Crown-of-thorns starfish */
# endif

# define AIR_SEA_GAS_EXCHANGE /*Original CPP flag */

/*** Coral Polyp model options. ***/
# if defined CORAL_POLYP
#  define CORAL_ZOOXANTHELLAE   /*Original CPP flag */
#  if defined ORGANIC_MATTER
#   define CORAL_MUCUS           /*Original CPP flag */
#   define CORAL_INGESTION       /*Original CPP flag */
#  endif
#  define CORAL_SIZE_DYNAMICS   /*Original CPP flag */
#  if defined CARBON_ISOTOPE
#   define CORAL_CARBON_ISOTOPE  /*Original CPP flag */
/*#   define CORAL_NONE_CO2_EQ*/
#  endif
#  if defined NUTRIENTS
/*#   define CORAL_NUTRIENTS*/  /*Original CPP flag */
#  endif
/*#  define CORAL_BORON_ISOTOPE*/  /*Original CPP flag */
# endif

#endif


/*** Box model option ***/

#define TESTMODE

#if defined CORAL_POLYP
# define CORAL_TESTMODE
/*# define CORAL_DEBUG*/
#endif

/*----------------------------------------------------*/

