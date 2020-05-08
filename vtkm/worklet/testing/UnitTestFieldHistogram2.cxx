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

#include <fstream>
#include <iostream>
#include <cstdio>
#include <vector>
//
// Make a simple 2D, 1000 point dataset populated with stat distributions
//
vtkm::cont::DataSet ReadTestDataSet()
{
  vtkm::cont::DataSet dataSet;
  
  const int dimension = 3;
  const int xVerts = 500;
  const int yVerts = 500;
  const int zVerts = 100;
  const int nVerts = xVerts * yVerts * zVerts;

  const int xCells = xVerts - 1;
  const int yCells = yVerts - 1;
  const int zCells = zVerts - 1;
  const int nCells = xCells * yCells * zCells;
  
  vtkm::Float32 *data = (vtkm::Float32 *)malloc(nVerts * sizeof(vtkm::Float32));
  std::cout << "QQ" << '\n';
  std::ifstream fileIn("/home/max/pf01.bin", std::ios::binary);
  //std::ofstream fileOut("/home/max/out.bin",std::ios::out | std::ios::binary | std::ios::app);
  float f;
  int i=0;
  
  while (fileIn.read(reinterpret_cast<char*>(&f), sizeof(float))){
      data[i++]=f;
      //fileOut.write(reinterpret_cast<const char*>(&f),sizeof(float));
  }
  //fileOut.close();
  
  vtkm::cont::ArrayHandleUniformPointCoordinates coordinates(vtkm::Id3(xVerts, yVerts, zVerts));
  dataSet.AddCoordinateSystem(vtkm::cont::CoordinateSystem("coordinates", coordinates));

  // Set point scalars
  dataSet.AddField(vtkm::cont::make_Field(
    "p_data", vtkm::cont::Field::Association::POINTS, data, nVerts, vtkm::CopyFlag::On));

  // Set cell scalars
  dataSet.AddField(vtkm::cont::make_Field(
    "c_data", vtkm::cont::Field::Association::CELL_SET, data, nCells, vtkm::CopyFlag::On));

  vtkm::cont::CellSetStructured<dimension> cellSet;

  //Set regular structure
  cellSet.SetPointDimensions(vtkm::make_Vec(xVerts, yVerts, zVerts));
  dataSet.SetCellSet(cellSet);
  free(data);
  return dataSet;
}

//
// Print the histogram result and tally
//
void PrintHistogram2(vtkm::cont::ArrayHandle<vtkm::Id> bins,
                    vtkm::Id numberOfBlocks,
                    vtkm::Id numberOfBins,
                    vtkm::cont::ArrayHandle<vtkm::Float32> delta,
                    vtkm::cont::ArrayHandle<vtkm::Float32> blockMin)
{
  vtkm::cont::ArrayHandle<vtkm::Id>::PortalConstControl binPortal = bins.GetPortalConstControl();
  vtkm::cont::ArrayHandle<vtkm::Float32>::PortalConstControl deltaPortal = delta.GetPortalConstControl();
  vtkm::cont::ArrayHandle<vtkm::Float32>::PortalConstControl blockMinPortal = blockMin.GetPortalConstControl();

  vtkm::Id sum = 0;
  for (vtkm::Id j = 0; j < numberOfBlocks; j++)
  {
    std::cout << "Block[" << j << "]:" << std::endl;
    for (vtkm::Id i = 0; i < numberOfBins; i++)
    {
      vtkm::Float64 lo = blockMinPortal.Get(j) + (static_cast<vtkm::Float64>(i) * deltaPortal.Get(j));
      vtkm::Float64 hi = lo + deltaPortal.Get(j);
      sum += binPortal.Get(j*numberOfBins+i);
      std::cout << "    BIN[" << i << "] Range[" << lo << ", " << hi << "] = " << binPortal.Get(j*numberOfBins+i)
                << std::endl;
    }
    std::cout << std::endl;
  }
  VTKM_TEST_ASSERT(test_equal(sum, 25000000), "Histogram not full");
}

//
// Create a dataset with known point data and cell data (statistical distributions)
// Extract arrays of point and cell fields
// Create output structure to hold histogram bins
// Run FieldHistogram filter
//
void TestFieldHistogram2()
{
  vtkm::Id xBnum = 20;
  vtkm::Id yBnum = 20;
  vtkm::Id zBnum = 20;
  auto temp = vtkm::make_Vec(xBnum, yBnum, zBnum);
  
  // Create the output bin array
  vtkm::Id numberOfBins = 10;
  vtkm::Range range;
  //vtkm::Float32 delta;
  vtkm::cont::ArrayHandle<vtkm::Float32> delta;
  vtkm::cont::ArrayHandle<vtkm::Float32> blockMin;
  vtkm::cont::ArrayHandle<vtkm::Id> bins;
  bins.Allocate(numberOfBins);
  delta.Allocate(xBnum*yBnum*zBnum);
  blockMin.Allocate(xBnum*yBnum*zBnum);
  // Data attached is the poisson distribution
  vtkm::cont::DataSet ds = ReadTestDataSet();

  // Get point data
  vtkm::cont::ArrayHandle<vtkm::Float32> p_data;
  ds.GetField("p_data").GetData().CopyTo(p_data);

  vtkm::worklet::FieldHistogram2 histogram;
  // Run data
  histogram.Run(p_data, temp, numberOfBins, blockMin, delta, bins);
  
  
  // std::cout << "Normal distributed POINT data:" << std::endl;
  // PrintHistogram2(bins, xBnum*yBnum*zBnum, numberOfBins, delta,blockMin);
  
  
} // TestFieldHistogram

int UnitTestFieldHistogram2(int argc, char* argv[])
{
  return vtkm::cont::testing::Testing::Run(TestFieldHistogram2, argc, argv);
}
