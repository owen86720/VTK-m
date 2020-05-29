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
#include <cstdlib> 
#include <ctime>

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
  
  // vtkm::Float32 data[nVerts] = {
  //   1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,
  //   33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64
  // };
  

  vtkm::Float32 *data = (vtkm::Float32 *)malloc(nVerts * sizeof(vtkm::Float32));
  std::ifstream fileIn("/home/max/pf01.bin", std::ios::binary);
  float f;
  int i=0;
  while (fileIn.read(reinterpret_cast<char*>(&f), sizeof(float))){
      data[i++]=f;
  }

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
  //free(data);
  return dataSet;
}

void Reconstruct(vtkm::cont::ArrayHandle<vtkm::Id> bins,
                vtkm::Id numberOfBlocks,
                vtkm::Id numberOfBins,
                vtkm::cont::ArrayHandle<vtkm::Float32> delta,
                vtkm::cont::ArrayHandle<vtkm::Float32> blockMin,
                vtkm::Vec<vtkm::Id,3> blockInfo,
                vtkm::Vec<vtkm::Id,3> dataInfo)
{
  srand( time(NULL) );
  int i=0;
  vtkm::Id R[numberOfBins+1];
  vtkm::Id P[numberOfBins];
  vtkm::Id xbIndex,ybIndex,zbIndex,index;
  vtkm::Float32 *data = (vtkm::Float32 *)malloc(dataInfo[0]*dataInfo[1]*dataInfo[2] * sizeof(vtkm::Float32));
  vtkm::cont::ArrayHandle<vtkm::Id>::PortalConstControl binPortal = bins.GetPortalConstControl();
  vtkm::cont::ArrayHandle<vtkm::Float32>::PortalConstControl deltaPortal = delta.GetPortalConstControl();
  vtkm::cont::ArrayHandle<vtkm::Float32>::PortalConstControl blockMinPortal = blockMin.GetPortalConstControl();
  // for(vtkm::Id i = 0; i < dataInfo[0]*dataInfo[1]*dataInfo[2]; i++)
  //     std::cout << data[i] << ' ';
  std::cout << std::endl;
  vtkm::Id sum = 0;
  for (vtkm::Id j = 0; j < numberOfBlocks; j++)
  {
    std::cout << "Block[" << j << "]:\n";
    zbIndex = j / (blockInfo[0]*blockInfo[1]);
    ybIndex = vtkm::FMod(static_cast<vtkm::Float64>(j) , static_cast<vtkm::Float64>(blockInfo[0]*blockInfo[1])) / blockInfo[0];
    xbIndex = vtkm::FMod(static_cast<vtkm::Float64>(j) , static_cast<vtkm::Float64>(blockInfo[0]*blockInfo[1])) - ybIndex * blockInfo[0];
    //std::cout << xbIndex << " " << ybIndex << " " << zbIndex << "\n";
    vtkm::Id offset =  static_cast<vtkm::Id>(zbIndex*dataInfo[0]*dataInfo[1]*dataInfo[2] / blockInfo[2]+ybIndex*dataInfo[0]*dataInfo[1] / blockInfo[1]+xbIndex*dataInfo[0] / blockInfo[0]);
    std::cout << "Offset: " << offset << "\n";

    R[0] = blockMinPortal.Get(j);
    for (vtkm::Id i = 0; i < numberOfBins; i++)
    {
      vtkm::Float64 lo = blockMinPortal.Get(j) + (static_cast<vtkm::Float64>(i) * deltaPortal.Get(j));
      vtkm::Float64 hi = lo + deltaPortal.Get(j);
      R[i+1] = hi;
      if(i==0)
        P[i] = binPortal.Get(j*numberOfBins+i);
      else{
        P[i] = P[i-1] + binPortal.Get(j*numberOfBins+i);
      }
      //sum += binPortal.Get(j*numberOfBins+i);
      std::cout << "    BIN[" << i << "] Range[" << lo << ", " << hi << "] = " << binPortal.Get(j*numberOfBins+i) << std::endl;
    }
    std::cout << std::endl;

    // for(vtkm::Id i = 0; i < numberOfBins; i++)
    //   std::cout << P[i] << ' ';
    // std::cout << std::endl;
    
    index = offset;
    for(vtkm::Id z = 0; z < dataInfo[2]/blockInfo[2]; z++){
      for(vtkm::Id y = 0; y < dataInfo[1]/blockInfo[1]; y++){
        for(vtkm::Id x = 0; x < dataInfo[0]/blockInfo[0]; x++){
          index = offset + z*dataInfo[0]*dataInfo[1] + y*dataInfo[0] + x;
          auto temp = rand() % P[numberOfBins-1];
          //std::cout << "temp: " << temp << " max: " << P[numberOfBins-1] <<"\n" ;
          //data[index] = index ;
          for (vtkm::Id i = 0; i < numberOfBins; i++){
            if(temp<P[i]){
              auto val = (R[i+1] - R[i]) * rand() / (RAND_MAX + 1.0) + R[i];
              data[index] = val;
              break;
            }
          }
        }
      }
    }
  }

  std::ofstream fileOut("/home/max/out3.bin",std::ios::out | std::ios::binary | std::ios::app);
  // for(vtkm::Id i = 0; i < dataInfo[0]*dataInfo[1]*dataInfo[2]; i++){
  //     std::cout << data[i] << ' ';
  //     fileOut.write(reinterpret_cast<const char*>(&data[i]),sizeof(float));
  //     if((i+1)%4==0)
  //       std::cout << std::endl;
  // }
  fileOut.close();
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
  //VTKM_TEST_ASSERT(test_equal(sum, 25000000), "Histogram not full");
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
  auto blockInfo = vtkm::make_Vec(xBnum, yBnum, zBnum);
  vtkm::Id xSize = 500;
  vtkm::Id ySize = 500;
  vtkm::Id zSize = 100;
  auto dataInfo = vtkm::make_Vec(xSize, ySize, zSize);
  
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
  histogram.Run(p_data, blockInfo, numberOfBins, blockMin, delta, bins);
  

  //PrintHistogram2(bins, xBnum*yBnum*zBnum, numberOfBins, delta,blockMin);
  Reconstruct(bins, xBnum*yBnum*zBnum, numberOfBins, delta,blockMin,blockInfo,dataInfo);
  
  
} // TestFieldHistogram

int UnitTestFieldHistogram2(int argc, char* argv[])
{
  return vtkm::cont::testing::Testing::Run(TestFieldHistogram2, argc, argv);
}
