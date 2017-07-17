#include <HE5_HdfEosDef.h>
#include <cfortHdf.h>

/* Given HE5 datatype, find out the corresponding HDF5 data type, and then
   use HDF5 intrinsic function H5Tget_size to figure out its size in byte,
   return this size */
int HE5Tget_size_ ( int fdatatype )
{ hid_t  cdatatype = FAIL;            /* return data type ID  */
  size_t size;
  cdatatype = HE5_EHconvdatatype( fdatatype );
  if( cdatatype == FAIL ) 
     size = -1;
  else 
    size = (int) H5Tget_size( cdatatype );
  return size;

} 

/* static float           fill_float32 = -0X1P+100;
   static double          fill_float64 = -0X1P+100; 
*/
static float           fill_float32 = -0X1P100;
static double          fill_float64 = -0X1P100;

int
r4Fill_( float *r4 )
{
  *r4 = fill_float32;
  return 0;
}

int
r8Fill_( double *r8 )
{
  *r8 = fill_float64;
  return 0;
}

FCALLSCFUN1( INT, HE5Tget_size_, HE5TGET_SIZE_, he5tget_size_, INT)
FCALLSCFUN1( INT, r4Fill_, R4FILL_, r4fill_, PFLOAT )
FCALLSCFUN1( INT, r8Fill_, R8FILL_, r8fill_, PDOUBLE )


