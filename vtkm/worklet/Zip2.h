//============================================================================
//  Copyright (c) Kitware, Inc.
//  All rights reserved.
//  See LICENSE.txt for details.
//
//  This software is distributed WITHOUT ANY WARRANTY; without even
//  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
//  PURPOSE.  See the above copyright notice for more information.
//============================================================================

#ifndef vtk_m_worklet_Zip2_h
#define vtk_m_worklet_Zip2_h

#include <vtkm/cont/ArrayHandle.h>
#include <vtkm/cont/ArrayHandleCounting.h>
#include <vtkm/worklet/DispatcherMapField.h>
#include <vtkm/worklet/WorkletMapField.h>

#include <vtkm/cont/Field.h>

namespace
{
// GCC creates false positive warnings for signed/unsigned char* operations.
// This occurs because the values are implicitly casted up to int's for the
// operation, and than  casted back down to char's when return.
// This causes a false positive warning, even when the values is within
// the value types range
#if defined(VTKM_GCC)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wconversion"
#endif // gcc

#if defined(VTKM_GCC)
#pragma GCC diagnostic pop
#endif // gcc
}

namespace vtkm
{
namespace worklet
{

//simple functor that prints basic statistics
class Zip2
{
public:
  // Calculate the adjacent difference between values in ArrayHandle
  class Unzipped : public vtkm::worklet::WorkletMapField
  {
  public:
    using ControlSignature = void(FieldIn zippedEle, FieldOut outputFirst, FieldOut outputSecond);
    using ExecutionSignature = void(_1, _2, _3);
    using InputDomain = _1;

    template <typename WholeArrayType, typename FieldType, typename FieldType2>
    VTKM_EXEC void operator()(const WholeArrayType& zipped,
                              FieldType& firstout,
                              FieldType2& secondout) const
    {
      firstout = zipped.first;
      secondout = zipped.second;
    }
  };

  template <typename FieldType, typename FieldType2, typename Storage>
  void Run(vtkm::cont::ArrayHandleZip<vtkm::cont::ArrayHandle<FieldType, Storage> , vtkm::cont::ArrayHandle<FieldType2, Storage> > zippedArray,
           vtkm::cont::ArrayHandle<FieldType>& Array1,
           vtkm::cont::ArrayHandle<FieldType2>& Array2)
  {
    // Difference between adjacent items is the bin count
    vtkm::worklet::DispatcherMapField<Unzipped> dispatcher;
    dispatcher.Invoke(zippedArray, Array1, Array2);

  }
};
}
} // namespace vtkm::worklet

#endif // vtk_m_worklet_FieldHistogram_h
