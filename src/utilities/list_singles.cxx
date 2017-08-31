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

/*!
  \file
  \ingroup utilities
  \author Kris Thielemans

  \brief Back project an image.

  \par Usage:
  \verbatim
  back_project output-filename image_to_back_project template_proj_data_file [backprojector-parfile ]\n"
  \endverbatim
  The template_proj_data_file will be used to get the scanner, mashing etc. details
  (its data will \e not be used, nor will it be overwritten).

  The default projector uses the ray-tracing matrix.
  \par Example parameter file for specifying the back projector
  \verbatim
  Back Projector parameters:=
    type := Matrix
      Back projector Using Matrix Parameters :=
        Matrix type := Ray Tracing
         Ray tracing matrix parameters :=
         End Ray tracing matrix parameters :=
        End Back Projector Using Matrix Parameters :=
  End:=
  \endverbatim
*/

#include "stir/data/SinglesRatesFromGEHDF5.h"
#include <iostream>

static void print_usage_and_exit()
{
  std::cerr<<"\nUsage:\nlist_singles _listmode_filename \n";

  exit(EXIT_FAILURE);
}


int 
main (int argc, char * argv[])
{
  using namespace stir;

  if (argc!=2)
    print_usage_and_exit();
  
  std::cout << "filename: " << argv[1] << std::endl;
  const std::string _listmode_filename = argv[1];

  SinglesRatesFromGEHDF5  singles;
  singles.read_singles_from_listmode_file(_listmode_filename);
  singles.write(std::cout);

  DetectionPosition<> pos(1,3,0);
  std::cout << "one of them: "<< singles.get_singles_rate(pos,1.,2.)<< std::endl;
  return EXIT_SUCCESS;
}

