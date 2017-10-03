/*
    Copyright (C) 2014, University College London
    This file is part of STIR.

    This file is free software; you can redistribute it and/or modify
    it under the terms of the Lesser GNU General Public License as published by
    the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.

    This file is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    Lesser GNU General Public License for more details.

    See STIR/LICENSE.txt for details

*/


#include "stir/data/SinglesRatesFromGEHDF5.h"
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
#include "stir/ProjData.h"
#include "stir/ProjDataInfoCylindricalNoArcCorr.h"
#include "stir/shared_ptr.h"
#include "stir/is_null_ptr.h"
#include "stir/Succeeded.h"
#include <numeric>
#include <fstream>
#include <iostream>

#ifndef STIR_NO_NAMESPACES
using std::cerr;
using std::endl;
using std
#include <iostream>

static void print_usage_and_exit()
{
  std::cerr<<"\nUsage:\nbin_to_GE _input_sinogram _output_sinogram\n";

  exit(EXIT_FAILURE);
}


int
main (int argc, char * argv[])
{
  using namespace stir;

  if (argc!=2)
    print_usage_and_exit();

  std::cout << "filename: " << argv[1] << std::endl;
  const std::string _input_sinogram = argv[1];

  const shared_ptr<ProjData> first_operand=ProjData::read_from_file(argv[1]);
  shared_ptr<ProjData> second_operand=first_operand->;
   shared_ptr<ProjDataInfo> pdi_ptr=first_operand->get_proj_data_info_sptr();
  Bin bin=pdi_ptr->get_bin();
  //PW need to chnage bin view number here and tangential position

  //Modify this bit as well.
  const std::string _output_sinogram = argv[2];
          shared_ptr<ProjData>
            out_proj_data_ptr(
                      new ProjDataInterfile(second_operand->get_exam_info_sptr(), second_operand->get_proj_data_info_ptr(), output_file_name, std::ios::in|std::ios::out|std::ios::trunc));



  return EXIT_SUCCESS;
}

