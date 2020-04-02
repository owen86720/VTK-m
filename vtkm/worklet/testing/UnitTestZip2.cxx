//============================================================================
//  Copyright (c) Kitware, Inc.
//  All rights reserved.
//  See LICENSE.txt for details.
//
//  This software is distributed WITHOUT ANY WARRANTY; without even
//  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
//  PURPOSE.  See the above copyright notice for more information.
//============================================================================

//#include <vtkm/worklet/FieldHistogram.h>

#include <vtkm/cont/ArrayHandle.h>
#include <vtkm/cont/Algorithm.h>
#include <vtkm/cont/testing/Testing.h>
#include <vtkm/worklet/Zip2.h>

void TestZip2()
{
  std::vector<vtkm::Id> inpBuffer{2,1,0,3,5,4,6,7,10,9};
  std::vector<vtkm::Id> keyBuffer{1,1,1,4,4,4,5,8,8,1};
  
  vtkm::cont::ArrayHandle<vtkm::Id> keys = vtkm::cont::make_ArrayHandle(keyBuffer);
  vtkm::cont::ArrayHandle<vtkm::Id> inps = vtkm::cont::make_ArrayHandle(inpBuffer);
  
  std::cout <<  "Before keys:";
  for (vtkm::Id i = 0; i < keys.GetNumberOfValues(); i++)
  {
    std::cout <<  keys.GetPortalConstControl().Get(i) << " ";
  }
  std::cout << std::endl << "Before inps:";
  for (vtkm::Id i = 0; i < inps.GetNumberOfValues(); i++)
  {
    std::cout <<  inps.GetPortalConstControl().Get(i) << " ";
  }

  vtkm::cont::ArrayHandleZip<vtkm::cont::ArrayHandle<vtkm::Id> , vtkm::cont::ArrayHandle<vtkm::Id> > zips = vtkm::cont::make_ArrayHandleZip(keys , inps);
  
  vtkm::cont::Algorithm::Sort(zips);

  vtkm::worklet::Zip2 unzip;
  unzip.Run(zips,keys,inps);

  std::cout << std::endl << "After keys:";
  for (vtkm::Id i = 0; i < keys.GetNumberOfValues(); i++)
  {
    std::cout <<  keys.GetPortalConstControl().Get(i) << " ";
  }
  std::cout << std::endl << "After inps:";
  for (vtkm::Id i = 0; i < inps.GetNumberOfValues(); i++)
  {
    std::cout <<  inps.GetPortalConstControl().Get(i) << " ";
  }
  std::cout << std::endl;

  std::cout <<"XXXXXX";
} 

int UnitTestZip2(int argc, char* argv[])
{
  return vtkm::cont::testing::Testing::Run(TestZip2, argc, argv);
}
