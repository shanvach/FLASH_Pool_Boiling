#include "io_ncmpi_name_to_id.h"

void FTOC(io_ncmpi_name_to_id)(const int * const pFileID,
			       const char name[],
			       const int * const pLen,
			       int *pVarID)
{
  const int fileID = *pFileID;
  const int len = *pLen;
  int err;
  char null_term_name[MAX_STRING_LENGTH+1];

  assert(len > 0 && len <= MAX_STRING_LENGTH);
  strncpy(null_term_name, name, len);
  null_term_name[len] = '\0';

  err = ncmpi_inq_varid(fileID, null_term_name, pVarID);
  assert(err == NC_NOERR);
}
