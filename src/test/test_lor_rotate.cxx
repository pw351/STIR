/*
    Copyright (C) 2016, UCL
    Copyright (C) 2016, University of Hull
    This file is part of STIR.

    This file is free software; you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.

    This file is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.
    See STIR/LICENSE.txt for details
*/
/*!
  \ingroup test
  \brief Test class for Time-Of-Flight
  \author Nikos Efthimiou
*/
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
#include "stir/ViewSegmentNumbers.h"
#include "stir/RelatedViewgrams.h"
//#include "stir/geometry/line_distances.h"
#include "stir/Succeeded.h"
#include "stir/shared_ptr.h"
#include "stir/RunTests.h"
#include "stir/Scanner.h"
#include "boost/lexical_cast.hpp"
#include "stir/Bin.h"
#include "stir/LORCoordinates.h"
#include "stir/round.h"
#include "stir/num_threads.h"
#include "stir/info.h"
#include "stir/warning.h"

START_NAMESPACE_STIR

//! A class to test the rotation of the lor in STIR
//! a key tight.
//! \author Palak Wadhwa
//!

class lor_rotate_Tests : public RunTests
{
public:
    void run_tests();

private:

    void test_lor_info();

    shared_ptr<Scanner> test_scanner_sptr;
    shared_ptr<ProjDataInfo> test_proj_data_info_sptr;
    shared_ptr<DiscretisedDensity<3, float> > test_discretised_density_sptr;
    shared_ptr<ProjMatrixByBin> test_proj_matrix_sptr;
    shared_ptr<ProjectorByBinPair> projector_pair_sptr;
    shared_ptr<DataSymmetriesForViewSegmentNumbers> symmetries_used_sptr;
};

void
lor_rotate_Tests::run_tests()
{
    // New Scanner
    test_scanner_sptr.reset(new Scanner(Scanner::PETMR_Signa));

    // New Proj_Data_Info
    test_proj_data_info_sptr.reset(ProjDataInfo::ProjDataInfoCTI(test_scanner_sptr,
                                                                 1,test_scanner_sptr->get_num_rings() -1,
                                                                 test_scanner_sptr->get_num_detectors_per_ring()/2,
                                                                 test_scanner_sptr->get_max_num_non_arccorrected_bins(),
                                                                 /* arc_correction*/false));

    test_lor_info();
//    test_tof_geometry_1();

    // New Discretised Density
    test_discretised_density_sptr.reset( new VoxelsOnCartesianGrid<float> (*test_proj_data_info_sptr, 1.f,
                                                                           CartesianCoordinate3D<float>(0.f, 0.f, 0.f),
                                                                           CartesianCoordinate3D<int>(-1, -1, -1)));
    // New ProjMatrix
    test_proj_matrix_sptr.reset(new ProjMatrixByBinUsingRayTracing());
    dynamic_cast<ProjMatrixByBinUsingRayTracing*>(test_proj_matrix_sptr.get())->set_num_tangential_LORs(1);
    dynamic_cast<ProjMatrixByBinUsingRayTracing*>(test_proj_matrix_sptr.get())->set_up(test_proj_data_info_sptr, test_discretised_density_sptr);
    shared_ptr<ForwardProjectorByBin> forward_projector_ptr(
                new ForwardProjectorByBinUsingProjMatrixByBin(test_proj_matrix_sptr));
    shared_ptr<BackProjectorByBin> back_projector_ptr(
                new BackProjectorByBinUsingProjMatrixByBin(test_proj_matrix_sptr));

    projector_pair_sptr.reset(
                new ProjectorByBinPairUsingSeparateProjectors(forward_projector_ptr, back_projector_ptr));
    projector_pair_sptr->set_up(test_proj_data_info_sptr, test_discretised_density_sptr);

    symmetries_used_sptr.reset(projector_pair_sptr->get_symmetries_used()->clone());

    // Deactivated it now because it takes a long time to finish.
    //        test_cache();
}

void
TOF_Tests::test_lor_info()
{

    ProjDataInfoCylindrical* proj_data_ptr =
            dynamic_cast<ProjDataInfoCylindrical*> (test_proj_data_info_sptr.get());
    Bin org_bin(0,15,0,24, /* value*/1);

           CartesianCoordinate3D<float> coord_1;
           CartesianCoordinate3D<float> coord_2;

           proj_data_ptr.find_cartesian_coordinates_of_detection(coord_1,coord_2,org_bin);

std::cout<<coord_1.x()<<"This is the x- coordinate of the detector" << std::endl;
std::cout<<coord_1.y()<<"This is the y- coordinate of the detector" << std::endl;
std::cout<<coord_1.z()<<"This is the z- coordinate of the detector" << std::endl;
std::cout<<coord_2.x()<<"This is the x- coordinate of the detector" << std::endl;
std::cout<<coord_2.y()<<"This is the y- coordinate of the detector" << std::endl;
std::cout<<coord_1.z()<<"This is the z- coordinate of the detector" << std::endl;


const int view2 = 15-test_scanner_sptr->get_default_intrinsic_tilt()/proj_data_ptr->get_azimuthal_angle_sampling();

Bin new_bin(0,view2,0,24,1);
CartesianCoordinate3D<float> coord_1_new;
CartesianCoordinate3D<float> coord_2_new;

    std::cerr<< lor_point_1.x() << " " << lor_point_1.y() << " " << lor_point_1.z() << " " <<
                lor_point_2.x() << " " << lor_point_2.y() << " " << lor_point_2.z() << std::endl;


    check_if_equal(correct_tof_mashing_factor,
                   test_proj_data_info_sptr->get_tof_mash_factor(), "Different TOF mashing factor.");

    check_if_equal(num_timing_positions,
                   test_proj_data_info_sptr->get_num_tof_poss(), "Different number of timing positions.");

    for (int timing_num = test_proj_data_info_sptr->get_min_tof_pos_num(), counter = 0;
         timing_num <= test_proj_data_info_sptr->get_max_tof_pos_num(); ++ timing_num, counter++)
    {
        Bin bin(0, 0, 0, 0, timing_num, 1.f);

        check_if_equal(static_cast<double>(correct_width_of_tof_bin),
                       static_cast<double>(test_proj_data_info_sptr->get_sampling_in_k(bin)), "Error in get_sampling_in_k()");
        check_if_equal(static_cast<double>(correct_timing_locations[counter]),
                       static_cast<double>(test_proj_data_info_sptr->get_k(bin)), "Error in get_sampling_in_k()");
    }

    float total_width = test_proj_data_info_sptr->get_k(Bin(0,0,0,0,test_proj_data_info_sptr->get_max_tof_pos_num(),1.f))
            - test_proj_data_info_sptr->get_k(Bin(0,0,0,0,test_proj_data_info_sptr->get_min_tof_pos_num(),1.f))
            + test_proj_data_info_sptr->get_sampling_in_k(Bin(0,0,0,0,0,1.f));

    set_tolerance(static_cast<double>(0.005));
    check_if_equal(static_cast<double>(total_width), static_cast<double>(test_proj_data_info_sptr->get_coincidence_window_width()),
                   "Coincidence widths don't match.");


}

END_NAMESPACE_STIR

int main()
{
    USING_NAMESPACE_STIR
    TOF_Tests tests;
    tests.run_tests();
    return tests.main_return_value();
}
