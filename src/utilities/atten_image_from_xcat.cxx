/*
 Copyright (C) 2011 - 2013, King's College London
 This file is part of STIR.

 This file is free software; you can redistribute it and/or modify
 it under the terms of the GNU Lesser General Public License as published by
 the Free Software Foundation; either version 2.3 of the License, or
 (at your option) any later version.

 This file is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU Lesser General Public License for more details.

 See STIR/LICENSE.txt for details
 */
/*!
 \file
 \ingroup utilities

 \brief This program Converts 1/pixel values of xcat phantom to Hounsfield units in cm-1 in order to get attenuation image.
 \author Palak Wadhwa
 */
#include "stir/DiscretisedDensity.h"
#include "stir/Succeeded.h"
#include "stir/IO/OutputFileFormat.h"
USING_NAMESPACE_STIR

int main(int argc, char **argv)
{
  if(argc!=3) {
    std::cerr << "Usage: " << argv[0]
              << "<output filename> <input filename>\n";
    exit(EXIT_FAILURE);
  }
  char const * const output_filename_prefix = argv[1];
  char const * const input_filename= argv[2];

  const std::auto_ptr<DiscretisedDensity<3,float> > image_aptr(DiscretisedDensity<3,float>::read_from_file(input_filename));
  const std::auto_ptr<DiscretisedDensity<3,float> > out_image_aptr(image_aptr->clone());

  DiscretisedDensity<3,float>& image = *image_aptr.get ();
  DiscretisedDensity<3,float>& output = *out_image_aptr.get ();

  const int min_z = image.get_min_index();
  const int max_z = image.get_max_index();


      for (int z=min_z; z<=max_z; z++)
        {

          const int min_y = image[z].get_min_index();
          const int max_y = image[z].get_max_index();



            for (int y=min_y;y<= max_y;y++)
              {

                const int min_x = image[z][y].get_min_index();
                const int max_x = image[z][y].get_max_index();



                  for (int x=min_x;x<= max_x;x++)
                  {

                             output[z][y][x]=4.97031154518263*image[z][y][x]-0.000170930400849;
                    }


                  }
              }

  //std::cout<<" image "<< image[0][-172][172]<<std::endl;
  const Succeeded res = OutputFileFormat<DiscretisedDensity<3,float> >::default_sptr()->
    write_to_file(output_filename_prefix, output);

  return res==Succeeded::yes ? EXIT_SUCCESS : EXIT_FAILURE;
}
