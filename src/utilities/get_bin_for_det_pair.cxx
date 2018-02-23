#include "stir/ProjDataInterfile.h"
#include "stir/DiscretisedDensity.h"
#include "stir/IO/read_from_file.h"
#include "stir/ProjDataFromStream.h"
#include "stir/SegmentByView.h"
#include "stir/SegmentBySinogram.h"
#include "stir/Sinogram.h"
#include "stir/Viewgram.h"

//#include "stir/Scanner.h"
#include "stir/ArrayFunction.h"
#include "stir/recon_array_functions.h"
#include "stir/display.h"
#include "stir/IO/interfile.h"
#include "stir/utilities.h"
#include "stir/shared_ptr.h"
#include "stir/is_null_ptr.h"

#include <numeric>
#include <fstream>
#include <iostream>
#include "stir/ProjDataInfoCylindricalNoArcCorr.h"
#include "stir/recon_buildblock/ProjMatrixByBin.h"
#include "stir/recon_buildblock/ProjMatrixByBinUsingRayTracing.h"
#include "stir/recon_buildblock/ForwardProjectorByBinUsingProjMatrixByBin.h"
#include "stir/recon_buildblock/BackProjectorByBinUsingProjMatrixByBin.h"
#include "stir/recon_buildblock/ProjMatrixElemsForOneBin.h"
#include "stir/recon_buildblock/ProjectorByBinPair.h"
#include "stir/recon_buildblock/ProjectorByBinPairUsingSeparateProjectors.h"
#include "stir/HighResWallClockTimer.h"
#include "stir/DiscretisedDensity.h"
#include "stir/VoxelsOnCartesianGrid.h"
#include "stir/recon_buildblock/ProjMatrixElemsForOneBin.h"
#include "stir/Succeeded.h"
#include <iostream>

static void print_usage_and_exit()
{
  std::cerr<<"\nUsage:\nget_bin_for_det_pair _projdata_filename \n";

  exit(EXIT_FAILURE);
}


int
main (int argc, char * argv[])
{
  using namespace stir;

  if (argc!=2)
    print_usage_and_exit();

  std::cout << "filename: " << argv[1] << std::endl;
  const std::string _projdata_filename = argv[1];

  shared_ptr<ProjData> proj_data_ptr =
    ProjData::read_from_file(argv[1]);

SegmentByView<float> seg0=proj_data_ptr->get_segment_by_view(0);
//shared_ptr<ProjDataInfoCylindricalNoArcCorr> pdi_ptr=seg0.get_proj_data_info_sptr();

Bin bin(0,0,0,0,0.f);
//Bin bin_pos=pdi_ptr->get_bin_for_det_pair(Bin,0,0,224,0);



}

