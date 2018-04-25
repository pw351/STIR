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
lor_rotate_Tests::test_lor_info()
{

    ProjDataInfoCylindricalNoArcCorr* proj_data_ptr =
            dynamic_cast<ProjDataInfoCylindricalNoArcCorr*> (test_proj_data_info_sptr.get());
    Bin org_bin(0,15,0,24, /* value*/1);

    int origDet1, origDet2;
    int origRing1, origRing2;

    int newDet1, newDet2;
    int newRing1, newRing2;

    proj_data_ptr->get_det_pair_for_bin(origDet1, origRing1, origDet2, origRing2,
                                        org_bin);
    // Test relationship between angle and view, through bins, for away of the center of FOV.
    for (int i = 0; i < proj_data_ptr->get_num_views() ; ++i)
    {
        Bin new_bin = org_bin;
        new_bin.view_num() += i;

        proj_data_ptr->get_det_pair_for_bin(newDet1, newRing1, newDet2, newRing2,
                                            new_bin);

        std::cout << /*"ODet1 " << origDet1 << " ODet2 " << origDet2 <<*/ "\n" <<
                     "NDet1 " << newDet1 << " NDet2 " << newDet2 << std::endl;

    }

    int nikos = 0;

//    CartesianCoordinate3D<float> coord_1;
//    CartesianCoordinate3D<float> coord_2;

//    proj_data_ptr->find_cartesian_coordinates_of_detection(coord_1,coord_2,org_bin);

//    int det1, det2,ring1, ring2;
//    proj_data_ptr->find_scanner_coordinates_given_cartesian_coordinates(det1,det2, ring1, ring2,coord_1,coord_2);
//    std::cout<<det1<<" This is the number of det1" << std::endl;
//    std::cout<<det2<<" This is the number of det2" << std::endl;


//    CartesianCoordinate3D<float> coord_1_new;
//    CartesianCoordinate3D<float> coord_2_new;

//    proj_data_ptr->find_cartesian_coordinates_of_detection(coord_1_new,coord_2_new,new_bin);

//    int det1_new, det2_new,ring1_new, ring2_new;
//    proj_data_ptr->find_scanner_coordinates_given_cartesian_coordinates(det1_new,det2_new, ring1_new, ring2_new,coord_1_new,coord_2_new);

//    float phi = test_scanner_sptr->get_default_intrinsic_tilt();
//    std::cout<<"The default intrinsic tilt is "<<phi<<std::endl;

//    CartesianCoordinate3D<float> coord_1_new_rotated;
//    coord_1_new_rotated.x() = coord_1.x()*cos(phi)-coord_1.y()*sin(phi);
//    coord_1_new_rotated.y() = coord_1.x()*sin(phi)+coord_1.y()*cos(phi);
//    coord_1_new_rotated.z() = coord_1.z();

//    CartesianCoordinate3D<float> coord_2_new_rotated;
//    coord_2_new_rotated.x() = coord_2.x()*cos(phi)-coord_2.y()*sin(phi);
//    coord_2_new_rotated.y() = coord_2.x()*sin(phi)+coord_2.y()*cos(phi);
//    coord_2_new_rotated.z() = coord_2.z();

//    std::cout<<coord_1_new_rotated.x()<<"This is the x- coordinate of the coord1_new_rotated" << std::endl;
//    std::cout<<coord_1_new_rotated.y()<<"This is the y- coordinate of the coord1_new_rotated" << std::endl;
//    std::cout<<coord_1_new_rotated.z()<<"This is the z- coordinate of the coord1_new_rotated" << std::endl;
//    std::cout<<coord_2_new_rotated.x()<<"This is the x- coordinate of the coord2_new_rotated" << std::endl;
//    std::cout<<coord_2_new_rotated.y()<<"This is the y- coordinate of the coord2_new_rotated" << std::endl;
//    std::cout<<coord_2_new_rotated.z()<<"This is the z- coordinate of the coord2_new_rotated" << std::endl;

}

END_NAMESPACE_STIR

int main()
{
    USING_NAMESPACE_STIR
            lor_rotate_Tests tests;
    tests.run_tests();
    return tests.main_return_value();
}
