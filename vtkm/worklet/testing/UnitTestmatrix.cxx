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

void Testmatrix()
{
  std::vector<vtkm::Id> inpBuffer{2,1,0,3,5,4,6,7,10,9};
  //inpBuffer = inpBuffer -1;
  std::vector<vtkm::Id> keyBuffer{1,1,1,4,4,4,5,8,8,1};
  
  vtkm::cont::ArrayHandle<vtkm::Id> keys = vtkm::cont::make_ArrayHandle(keyBuffer);
  vtkm::cont::ArrayHandle<vtkm::Id> inps = vtkm::cont::make_ArrayHandle(inpBuffer);
  
  for (vtkm::Id i = 0; i < keys.GetNumberOfValues(); i++)
  {
    std::cout <<  keys.GetPortalConstControl().Get(i) << " ";
  }
  for (vtkm::Id i = 0; i < inps.GetNumberOfValues(); i++)
  {
    std::cout <<  inps.GetPortalConstControl().Get(i) << " ";
  }

  vtkm::cont::ArrayHandleZip<vtkm::cont::ArrayHandle<vtkm::Id> , vtkm::cont::ArrayHandle<vtkm::Id> > zips = vtkm::cont::make_ArrayHandleZip(keys , inps);

  vtkm::cont::ArrayHandle<vtkm::Id> uniqueKeys;
  vtkm::cont::ArrayHandle<vtkm::Id> sums;
  vtkm::cont::ArrayHandle<vtkm::Id> products;

  //vtkm::cont::Algorithm::ReduceByKey(keys,inps,uniqueKeys,sums,vtkm::Subtract());
  //vtkm::cont::Algorithm::ReduceByKey(keys,inps,uniqueKeys,products,vtkm::Multiply());
  
  //vtkm::cont::Algorithm::Sort(zips);


  std::cout << "XXXXX" << std::endl;

  vtkm::worklet::Zip2 unzip;
  unzip.Run(zips,sums,products);

  //vtkm::cont::Algorithm::SortByKey(keys,inps);
  /*
  for (vtkm::Id i = 0; i < keys.GetNumberOfValues(); i++)
  {
    std::cout <<  keys.GetPortalConstControl().Get(i) << " ";
  }
  for (vtkm::Id i = 0; i < inps.GetNumberOfValues(); i++)
  {
    std::cout <<  inps.GetPortalConstControl().Get(i) << " ";
  }
  */
  for (vtkm::Id i = 0; i < sums.GetNumberOfValues(); i++)
  {
    std::cout <<  sums.GetPortalConstControl().Get(i) << " ";
  }
  for (vtkm::Id i = 0; i < products.GetNumberOfValues(); i++)
  {
    std::cout <<  products.GetPortalConstControl().Get(i) << " ";
  }
  /*
  vtkm::cont::ArrayHandle<vtkm::Id>::PortalConstControl SumPortal = sums.GetPortalConstControl();
  vtkm::cont::ArrayHandle<vtkm::Id>::PortalConstControl ProductPortal = products.GetPortalConstControl();
  
  std::cout << "SUM:" << sums.GetNumberOfValues() << std::endl; 
  for (vtkm::Id i = 0; i < sums.GetNumberOfValues(); i++)
  {
    std::cout << "sum[" << i << "]:" <<  SumPortal.Get(i) << std::endl ;
  }
  
  std::cout << "PRODUCT:" << products.GetNumberOfValues() << std::endl;
  for (vtkm::Id i = 0; i < products.GetNumberOfValues(); i++)
  {
    std::cout << "product[" << i << "]:" <<  ProductPortal.Get(i) << std::endl ;
  }
  */
  std::cout << "FQTTTTTTTTTTTt:" << std::endl;
} 

int UnitTestmatrix(int argc, char* argv[])
{
  return vtkm::cont::testing::Testing::Run(Testmatrix, argc, argv);
}
