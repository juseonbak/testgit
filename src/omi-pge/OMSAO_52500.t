# Status Message Text for Toolkit --- OMI SAO PGE
#
# Error messages for input module
#
# Author:   T.P. Kurosu
# Date:     31 December, 2002
# Version:  0.5
# History: 
#
%INSTR = OMI
%LABEL = OMSAO
%SEED  = 52500

OMSAO_W_GETLUN    WARNING...failed to read LUN from PCF
OMSAO_E_GETLUN    ERROR...failed to read LUN from PCF. PGE aborting with exit code 1

OMSAO_E_GET_MOLINDEX      ERROR...failed to convert molecule ID to index. PGE aborting with exit code 1
OMSAO_S_GET_MOLINDEX      SUCCESS...identified molecule ID for current PGE run

OMSAO_E_OPEN_FITCTRL_FILE  ERROR...failed to open fitting control file. PGE aborting with exit code 1
OMSAO_E_READ_FITCTRL_FILE  ERROR...failed during READ from fitting control file. PGE aborting with exit code 1
OMSAO_W_READ_FITCTRL_FILE  WARNING...failed during READ from fitting control file.
OMSAO_W_CLOSE_FITCTRL_FILE WARNING...failed to close fitting control file
OMSAO_S_READ_FITCTRL_FILE  SUCCESS...completed READ from fitting control file

OMSAO_E_GET_MOLFITNAME     ERROR...failed to convert mol-name to fit-index. PGE aborting with exit code 1

OMSAO_E_OPEN_REFSPEC_FILE  ERROR...failed to open reference spectrum file. PGE aborting with exit code 1
OMSAO_E_READ_REFSPEC_FILE  ERROR...failed to READ from reference spectrum file. PGE aborting with exit code 1
OMSAO_W_CLOSE_REFSPEC_FILE WARNING...failed to close reference spectrum file
OMSAO_S_READ_REFSPEC_FILE  SUCCESS...completed READ of all reference spectra

OMSAO_W_READ_AMFTABLE_NAME  WARNING...failed to read name of AMF table file from PCF
OMSAO_W_OPEN_AMFTABLE_FILE  WARNING...failed to open AMF table file
OMSAO_W_READ_AMFTABLE_FILE  WARNING...failed during READ from AMF table file
OMSAO_W_CLOSE_AMFTABLE_FILE WARNING...failed to close AMF table file
OMSAO_S_READ_AMFTABLE_FILE  SUCCESS...completed READ of AMF table file
OMSAO_W_AMFTABLE_MEMALLOC   WARNING...failed to allocate memory for AMF table

OMSAO_E_OPEN_RANCIL_FILE   ERROR...failed to open ancillary file for reading. PGE aborting with exit code 1
OMSAO_S_OPEN_RANCIL_FILE   SUCCESS...opened ancillary file for reading
OMSAO_E_OPEN_WANCIL_FILE   ERROR...failed to open ancillary file for writing. PGE aborting with exit code 1
OMSAO_S_OPEN_WANCIL_FILE   SUCCESS...openened ancillary file for writing
OMSAO_E_OPEN_UANCIL_FILE   ERROR...unknown action in OPEN of ancillary file
OMSAO_E_UACTION_ANCIL_FILE ERROR...unknown request to anciallry file
OMSAO_W_CLOSE_ANCIL_FILE   WARNING...failed to close ancillary file
OMSAO_S_CLOSE_ANCIL_FILE   SUCCESS...closed ancillary file

OMSAO_E_OPEN_L1B_FILE      ERROR...failed to open L1b file. PGE aborting with exit code 1
OMSAO_E_READ_L1B_FILE      ERROR...failed to read from L1b file. PGE aborting with exit code 1
OMSAO_W_CLOS_L1B_FILE      WARNING...failed to close OMI L1B file

OMSAO_E_OPEN_L2CLD_FILE      ERROR...failed to open OMICLD L2 file. PGE aborting with exit code 1
OMSAO_E_READ_L2CLD_FILE      ERROR...failed to read from OMICLD L2 file. PGE aborting with exit code 1
OMSAO_W_CLOS_L2CLD_FILE      WARNING...failed to close OMICLD L2 file

OMSAO_E_MEM_ALLOC          ERROR...memory allocation failure. PGE aborting with exit code 1
OMSAO_E_MEM_DALLOC         ERROR...memory de-allocation failure. PGE aborting with exit code 1

OMSAO_E_INTERPOL           ERROR...interpolation failed. PGE aborting with exit code 1
OMSAO_W_INTERPOL           WARNING...interpolation failed
OMSAO_E_INTERPOL_REFSPEC   ERROR...interpolation failed for reference spectum. PGE aborting with exit code 1
OMSAO_W_INTERPOL_REFSPEC   WARNING...interpolation found incomplete reference spectrum for:

OMSAO_E_HE5SWOPEN          ERROR...failed to open HE5 output file. PGE aborting with exit code 1
OMSAO_E_HE5SWCREATE        ERROR...failed to create HE5 swath. PGE aborting with exit code 1
OMSAO_E_HE5SWDEFDIM        ERROR...failed to define HE5 swath dimensions. PGE aborting with exit code 1
OMSAO_E_HE5SWDEFFLD        ERROR...failed to define HE5 swath geo/data fields. PGE aborting with exit code 1
OMSAO_E_HE5SWATTACH        ERROR...failed to attach to existing HE5 swath. PGE aborting with exit code 1
OMSAO_E_HE5SWWRFLD         ERROR...failed to write to HE5 data field. PGE aborting with exit code 1
OMSAO_W_HE5SWCLOSE         WARNING...failed to close HE5 output file

OMSAO_E_GETATTR            ERROR...failed to get MetaData Attribute. PGE aborting with exit code 1
OMSAO_E_METINIT            ERROR...failed to initialize MCF. PGE aborting with exit code 1
OMSAO_E_GETREF             ERROR...failed to get L2 Input Pointer. PGE aborting with exit code 1
OMSAO_E_MDMISMATCH         ERROR...mismatch for Metadata value field. PGE aborting with exit code 1
OMSAO_E_LGIDCOMP           ERROR...could not compose LocalGranuleID. PGE aborting with exit code 1
OMSAO_E_MDL2INV            ERROR...failed to set L2 inventory Metadata. PGE aborting with exit code 1
OMSAO_E_MDINIT             ERROR...failed to initiate L2 Metadata. PGE aborting with exit code 1
OMSAO_E_CMDWRT             ERROR...failed to write CoreMetaData. PGE aborting with exit code 1
OMSAO_W_MDMISMATCH         WARNING...mismatch for Metadata value field. PGE continuing with exit code 1
OMSAO_W_MDL2INV            WARNING...failed to set L2 inventory Metadata.
OMSAO_W_MDL2ARC            WARNING...failed to set L2 archived Metadata.
OMSAO_W_GETATTR            WARNING...failed to get MetaData Attribute.

OMSAO_S_ENDOFRUN           END_OF_RUN_S...PGE execution completed normally, exit code 0
OMSAO_W_ENDOFRUN           END_OF_RUN_W...PGE execution completed with warnings, exit code 1
OMSAO_E_ENDOFRUN           END_OF_RUN_E...PGE execution aborted abnormally, exit code 1
OMSAO_F_ENDOFRUN           END_OF_RUN_F...PGE execution aborted abnormally, exit code 1
OMSAO_U_ENDOFRUN           END_OF_RUN_U...PGE execution completed with unknown exit code
