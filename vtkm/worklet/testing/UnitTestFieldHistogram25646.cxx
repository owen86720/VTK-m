//============================================================================
//  Copyright (c) Kitware, Inc.
//  All rights reserved.
//  See LICENSE.txt for details.
//
//  This software is distributed WITHOUT ANY WARRANTY; without even
//  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
//  PURPOSE.  See the above copyright notice for more information.
//============================================================================

#include <vtkm/worklet/FieldHistogram2.h>

#include <vtkm/cont/DataSet.h>
#include <vtkm/cont/testing/Testing.h>

//
// Make a simple 2D, 1000 point dataset populated with stat distributions
//
vtkm::cont::DataSet MakeTestDataSet2()
{
  vtkm::cont::DataSet dataSet;

  const int dimension = 2;
  const int xVerts = 20;
  const int yVerts = 50;
  const int nVerts = xVerts * yVerts;

  const int xCells = xVerts - 1;
  const int yCells = yVerts - 1;
  const int nCells = xCells * yCells;

  // Normal distribution [0:49] mean = 25 standard deviation = 5.0
  vtkm::Float32 normal[nVerts] = {
    24, 19, 28, 19, 25, 28, 25, 22, 27, 26, 35, 26, 30, 28, 24, 23, 21, 31, 20, 11, 21, 22, 14, 25,
    20, 24, 24, 21, 24, 29, 26, 21, 32, 29, 23, 28, 31, 25, 23, 30, 18, 24, 22, 25, 33, 24, 22, 23,
    21, 17, 20, 28, 30, 18, 20, 32, 25, 24, 32, 15, 27, 24, 27, 19, 30, 27, 17, 24, 29, 23, 22, 19,
    24, 19, 28, 24, 25, 24, 25, 30, 24, 31, 30, 27, 25, 25, 25, 15, 29, 23, 29, 29, 21, 25, 35, 24,
    28, 10, 31, 23, 22, 22, 22, 33, 29, 27, 18, 27, 27, 24, 20, 20, 21, 29, 23, 31, 23, 23, 22, 23,
    30, 27, 28, 31, 16, 29, 25, 19, 33, 28, 25, 24, 15, 27, 37, 29, 15, 19, 14, 19, 24, 23, 30, 29,
    35, 22, 19, 26, 26, 14, 24, 30, 32, 23, 30, 29, 26, 27, 25, 23, 17, 26, 32, 29, 20, 17, 21, 23,
    22, 20, 36, 12, 26, 23, 15, 29, 24, 22, 26, 33, 24, 23, 20, 26, 22, 17, 26, 26, 34, 22, 26, 17,
    23, 18, 29, 27, 21, 29, 28, 29, 24, 25, 28, 19, 18, 21, 23, 23, 27, 25, 24, 25, 24, 25, 21, 25,
    21, 27, 23, 20, 29, 15, 28, 30, 24, 27, 17, 23, 16, 21, 25, 17, 27, 28, 21, 13, 19, 27, 16, 30,
    31, 25, 30, 17, 17, 25, 26, 22, 21, 17, 24, 17, 25, 22, 27, 14, 27, 24, 27, 25, 26, 31, 21, 23,
    30, 30, 22, 19, 23, 22, 23, 25, 24, 25, 24, 28, 26, 30, 18, 25, 30, 37, 27, 34, 28, 34, 25, 10,
    25, 22, 35, 30, 24, 32, 24, 34, 19, 29, 26, 16, 27, 17, 26, 23, 27, 25, 26, 21, 31, 21, 28, 15,
    32, 24, 23, 23, 18, 15, 22, 25, 16, 25, 31, 26, 25, 28, 24, 26, 23, 25, 33, 20, 27, 28, 24, 29,
    32, 20, 24, 20, 19, 32, 24, 6,  24, 21, 26, 18, 15, 30, 19, 26, 22, 30, 35, 23, 22, 30, 20, 22,
    18, 30, 28, 25, 16, 25, 27, 30, 18, 24, 30, 28, 20, 19, 20, 28, 21, 24, 15, 33, 20, 18, 20, 36,
    30, 26, 25, 18, 28, 27, 31, 31, 15, 26, 16, 22, 27, 14, 17, 27, 27, 22, 32, 30, 22, 34, 22, 25,
    20, 22, 26, 29, 28, 33, 18, 23, 20, 20, 27, 24, 28, 21, 25, 27, 25, 19, 19, 25, 19, 32, 29, 27,
    23, 21, 28, 33, 23, 23, 28, 26, 31, 19, 21, 29, 21, 27, 23, 32, 24, 26, 21, 28, 28, 24, 17, 31,
    27, 21, 19, 32, 28, 23, 30, 23, 29, 15, 26, 26, 15, 20, 25, 26, 27, 31, 21, 23, 23, 33, 28, 19,
    23, 22, 22, 25, 27, 17, 23, 17, 25, 28, 26, 30, 32, 31, 19, 25, 25, 19, 23, 29, 27, 23, 34, 22,
    13, 21, 32, 10, 20, 33, 21, 17, 29, 31, 14, 24, 23, 19, 19, 22, 17, 26, 37, 26, 22, 26, 38, 29,
    29, 27, 30, 20, 31, 14, 32, 32, 24, 23, 23, 18, 21, 31, 24, 20, 28, 15, 21, 25, 25, 20, 30, 25,
    22, 21, 21, 25, 24, 25, 18, 23, 28, 30, 20, 27, 27, 19, 10, 32, 24, 20, 29, 26, 25, 20, 25, 29,
    28, 24, 32, 26, 22, 19, 23, 27, 27, 29, 20, 25, 21, 30, 28, 31, 24, 19, 23, 19, 19, 18, 30, 18,
    16, 24, 20, 20, 30, 25, 29, 25, 31, 21, 28, 31, 24, 26, 27, 21, 24, 23, 26, 18, 32, 26, 28, 26,
    24, 26, 29, 30, 22, 20, 24, 28, 25, 29, 20, 21, 22, 15, 30, 27, 33, 26, 22, 32, 30, 31, 20, 19,
    24, 26, 27, 31, 17, 17, 33, 27, 16, 27, 27, 22, 27, 19, 24, 21, 17, 24, 28, 23, 26, 24, 19, 26,
    20, 24, 22, 19, 22, 21, 21, 28, 29, 39, 19, 16, 25, 29, 31, 22, 22, 29, 26, 22, 22, 22, 26, 23,
    23, 23, 30, 25, 25, 25, 27, 29, 18, 33, 21, 12, 22, 29, 12, 20, 35, 22, 34, 28, 18, 29, 21, 20,
    24, 33, 24, 26, 23, 34, 31, 25, 31, 22, 35, 21, 20, 29, 27, 22, 30, 22, 27, 23, 22, 32, 16, 19,
    27, 22, 24, 27, 21, 33, 25, 25, 19, 28, 20, 27, 21, 25, 28, 20, 27, 22, 21, 20, 26, 30, 33, 23,
    20, 24, 17, 23, 28, 35, 14, 23, 22, 28, 28, 26, 25, 18, 20, 28, 28, 22, 13, 24, 22, 20, 30, 26,
    26, 18, 22, 20, 23, 24, 20, 27, 34, 28, 18, 24, 34, 33, 25, 33, 37, 21, 20, 31, 19, 23, 29, 22,
    21, 24, 19, 27, 19, 32, 25, 23, 33, 26, 33, 27, 29, 30, 19, 22, 30, 19, 18, 24, 25, 17, 31, 19,
    31, 26, 22, 23, 28, 28, 25, 24, 19, 19, 27, 28, 23, 21, 29, 26, 31, 22, 22, 25, 16, 29, 21, 22,
    23, 25, 22, 21, 22, 19, 27, 26, 28, 30, 22, 21, 24, 22, 23, 26, 28, 22, 18, 25, 23, 27, 31, 19,
    15, 29, 20, 19, 27, 25, 21, 29, 22, 24, 25, 17, 36, 29, 22, 22, 24, 28, 27, 22, 26, 31, 29, 31,
    18, 25, 23, 16, 37, 27, 21, 31, 25, 24, 20, 23, 28, 33, 24, 21, 26, 20, 18, 31, 20, 24, 23, 19,
    27, 17, 23, 23, 20, 26, 28, 23, 26, 31, 25, 31, 19, 32, 26, 18, 19, 29, 20, 21, 15, 25, 27, 29,
    22, 22, 22, 26, 23, 22, 23, 29, 28, 20, 21, 22, 20, 22, 27, 25, 23, 32, 23, 20, 31, 20, 27, 26,
    34, 20, 22, 36, 21, 29, 25, 20, 21, 22, 29, 29, 25, 22, 24, 22
  };

  vtkm::cont::ArrayHandleUniformPointCoordinates coordinates(vtkm::Id3(xVerts, yVerts, 1));
  dataSet.AddCoordinateSystem(vtkm::cont::CoordinateSystem("coordinates", coordinates));

  // Set point scalars
  dataSet.AddField(vtkm::cont::make_Field(
    "p_normal", vtkm::cont::Field::Association::POINTS, normal, nVerts, vtkm::CopyFlag::On));

  // Set cell scalars
  dataSet.AddField(vtkm::cont::make_Field("c_normal",
                                          vtkm::cont::Field::Association::CELL_SET,
                                          "cells",
                                          normal,
                                          nCells,
                                          vtkm::CopyFlag::On));

  vtkm::cont::CellSetStructured<dimension> cellSet("cells");

  //Set regular structure
  cellSet.SetPointDimensions(vtkm::make_Vec(xVerts, yVerts));
  dataSet.AddCellSet(cellSet);

  return dataSet;
}

//
// Print the histogram result and tally
//
void PrintHistogram2(vtkm::cont::ArrayHandle<vtkm::Id> bins,
                    vtkm::Id numberOfBins,
                    const vtkm::Range& range,
                    vtkm::Float32 delta)
{
  vtkm::cont::ArrayHandle<vtkm::Id>::PortalConstControl binPortal = bins.GetPortalConstControl();

  vtkm::Id sum = 0;
  for (vtkm::Id i = 0; i < numberOfBins; i++)
  {
    vtkm::Float64 lo = range.Min + (static_cast<vtkm::Float64>(i) * delta);
    vtkm::Float64 hi = lo + delta;
    sum += binPortal.Get(i);
    std::cout << "  BIN[" << i << "] Range[" << lo << ", " << hi << "] = " << binPortal.Get(i)
              << std::endl;
  }
  VTKM_TEST_ASSERT(test_equal(sum, 1000), "Histogram not full");
}

//
// Create a dataset with known point data and cell data (statistical distributions)
// Extract arrays of point and cell fields
// Create output structure to hold histogram bins
// Run FieldHistogram filter
//
void TestFieldHistogram2()
{
  // Create the output bin array
  vtkm::Id numberOfBlocks = 3;
  vtkm::Id numberOfBins = 10;
  vtkm::Range range;
  vtkm::Float32 delta;
  vtkm::cont::ArrayHandle<vtkm::Id> bins;
  bins.Allocate(numberOfBins);

  // Data attached is the poisson distribution
  vtkm::cont::DataSet ds = MakeTestDataSet2();

  // Get point data
  vtkm::cont::ArrayHandle<vtkm::Float32> p_normal;
  ds.GetField("p_normal").GetData().CopyTo(p_normal);

  vtkm::worklet::FieldHistogram2 histogram;
  // Run data
  histogram.Run(p_normal, numberOfBlocks ,numberOfBins, range, delta, bins);
  std::cout << "Normal distributed POINT data:" << std::endl;
  PrintHistogram2(bins, numberOfBins, range, delta);

} // TestFieldHistogram

int UnitTestFieldHistogram2(int argc, char* argv[])
{
  return vtkm::cont::testing::Testing::Run(TestFieldHistogram2, argc, argv);
}
