//============================================================================
//  Copyright (c) Kitware, Inc.
//  All rights reserved.
//  See LICENSE.txt for details.
//
//  This software is distributed WITHOUT ANY WARRANTY; without even
//  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
//  PURPOSE.  See the above copyright notice for more information.
//============================================================================

#ifndef vtk_m_worklet_FieldHistogram2_h
#define vtk_m_worklet_FieldHistogram2_h

#include <vtkm/Math.h>
#include <vtkm/cont/Algorithm.h>
#include <vtkm/cont/ArrayGetValues.h>
#include <vtkm/cont/ArrayCopy.h>
#include <vtkm/cont/ArrayHandle.h>
#include <vtkm/cont/ArrayHandleCounting.h>
#include <vtkm/worklet/DispatcherMapField.h>
#include <vtkm/worklet/WorkletMapField.h>
#include <vtkm/worklet/Zip2.h>
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
template <typename T>
T compute_delta(T fieldMinValue, T fieldMaxValue, vtkm::Id num)
{
  using VecType = vtkm::VecTraits<T>;
  const T fieldRange = fieldMaxValue - fieldMinValue;
  return fieldRange / static_cast<typename VecType::ComponentType>(num);
}
#if defined(VTKM_GCC)
#pragma GCC diagnostic pop
#endif // gcc
}

namespace vtkm
{
namespace worklet
{

//simple functor that prints basic statistics
class FieldHistogram2
{
public:
  
  template <typename FieldType>
  class ComputeDelta : public vtkm::worklet::WorkletMapField
  {
  public:
    using ControlSignature = void(FieldIn blockMinValue, FieldIn blockMaxValue, FieldOut delta);
    using ExecutionSignature = void(_1, _2, _3);
    using InputDomain = _1;

    vtkm::Id NumberOfBins ;

    VTKM_CONT
    ComputeDelta(vtkm::Id NumberOfBins0)
      : NumberOfBins(NumberOfBins0)
    {
    }

    VTKM_EXEC
    void operator()(const FieldType& blockMinValue, const FieldType& blockMaxValue, FieldType& delta) const
    {
      using VecType = vtkm::VecTraits<FieldType>;
      const FieldType fieldRange = blockMaxValue - blockMinValue;
      delta = fieldRange / static_cast<typename VecType::ComponentType>(NumberOfBins);
    }
  };

  // For each value set the block it should be in
  class SetBlock : public vtkm::worklet::WorkletMapField
  {
  public:
    using ControlSignature = void(FieldIn index, FieldOut blockIndex);
    using ExecutionSignature = void(_1, _2);
    using InputDomain = _1;

    vtkm::Vec<vtkm::Id,3> cutInfo;
    vtkm::Id xVerts = 500 ;
    vtkm::Id yVerts = 500 ;
    vtkm::Id zVerts = 100 ;

    VTKM_CONT
    SetBlock(vtkm::Vec<vtkm::Id,3> cutInfo0)
      : cutInfo(cutInfo0)
    {
    }

    VTKM_EXEC
    void operator()(const vtkm::Id& index, vtkm::Id& blockIndex) const
    {
      vtkm::Id xIndex,yIndex,zIndex;
      vtkm::Id xbIndex,ybIndex,zbIndex;

      zIndex = index / (xVerts*yVerts);
      yIndex = vtkm::FMod(static_cast<vtkm::Float64>(index) , static_cast<vtkm::Float64>(xVerts*yVerts)) / xVerts;
      xIndex = vtkm::FMod(static_cast<vtkm::Float64>(index) , static_cast<vtkm::Float64>(xVerts*yVerts)) - yIndex * xVerts;

      xbIndex = xIndex/(xVerts/cutInfo[0]);
      ybIndex = yIndex/(yVerts/cutInfo[1]);
      zbIndex = zIndex/(zVerts/cutInfo[2]);

      blockIndex = xbIndex + ybIndex * cutInfo[0] + zbIndex * cutInfo[0] * cutInfo[1] ;

      // if(index == 24999999)
      //   std::cout << "FKKKK" << '\n';
      // if(index == 25000000){
      //   std::cout << "FKKKK222222" << '\n';
      //   std::cout << index << ' '<< zIndex <<' ' << yIndex <<' '<< xIndex <<' '<< zbIndex<<' ' << ybIndex <<' '<< xbIndex << '\n';
      // }
    }
  };

  // For each value set the bin it should be in
  template <typename FieldType>
  class SetHistogramBin : public vtkm::worklet::WorkletMapField
  {
  public:
    using ControlSignature = void(FieldIn blockIndex, FieldIn value, FieldOut binIndex);
    using ExecutionSignature = void(_1, _2, _3);
    using InputDomain = _1;

    vtkm::Id numberOfBins;
    vtkm::cont::ArrayHandle<vtkm::Float32> blockMin;
    vtkm::cont::ArrayHandle<vtkm::Float32> delta;
    VTKM_CONT
    SetHistogramBin(vtkm::Id numberOfBins0,
                    vtkm::cont::ArrayHandle<vtkm::Float32> blockMin0,
                    vtkm::cont::ArrayHandle<vtkm::Float32> delta0)
      : numberOfBins(numberOfBins0)
    {
      vtkm::cont::ArrayCopy(blockMin0 , blockMin );
      vtkm::cont::ArrayCopy(delta0 , delta );
    }

    VTKM_EXEC
    void operator()(const vtkm::Id& blockIndex, 
                    const FieldType& value, 
                    vtkm::Id& binIndex) const
    {
      FieldType minValue;
      FieldType deltaValue;
      minValue = blockMin.GetPortalConstControl().Get(blockIndex);
      deltaValue = delta.GetPortalConstControl().Get(blockIndex);

      binIndex = static_cast<vtkm::Id>((value - minValue) / deltaValue);
      if (binIndex < 0)
        binIndex = 0;
      else if (binIndex >= numberOfBins)
        binIndex = numberOfBins - 1;
    }
  };

  // Calculate the adjacent difference between values in ArrayHandle
  class BlockShift : public vtkm::worklet::WorkletMapField
  {
  public:
    using ControlSignature = void(FieldIn blockIndex, FieldIn binIndex, FieldOut output);
    using ExecutionSignature = void(_1, _2, _3);
    using InputDomain = _1;

    vtkm::Id numberOfBins;

    VTKM_CONT
    BlockShift(vtkm::Id numberOfBins0)
      : numberOfBins(numberOfBins0)
    {
    }
    
    VTKM_EXEC
    void operator()(const vtkm::Id& blockIndex,
                              const vtkm::Id& binIndex,
                              vtkm::Id& Added) const
    {
      Added = blockIndex*numberOfBins + binIndex;
    }
  };

  // Calculate the adjacent difference between values in ArrayHandle
  class AdjacentDifference : public vtkm::worklet::WorkletMapField
  {
  public:
    using ControlSignature = void(FieldIn inputIndex, WholeArrayIn counts, FieldOut outputCount);
    using ExecutionSignature = void(_1, _2, _3);
    using InputDomain = _1;

    template <typename WholeArrayType>
    VTKM_EXEC void operator()(const vtkm::Id& index,
                              const WholeArrayType& counts,
                              vtkm::Id& difference) const
    {
      if (index == 0)
        difference = counts.Get(index);
      else
        difference = counts.Get(index) - counts.Get(index - 1);
    }
  };

  // Execute the histogram binning filter given data and number of bins
  // Returns:
  // min value of the bins
  // delta/range of each bin
  // number of values in each bin
  // template <typename FieldType, typename Storage>
  // void Run(vtkm::cont::ArrayHandle<FieldType, Storage> fieldArray,
  //          vtkm::Vec<vtkm::Id,3> cutInfo,
  //          vtkm::Id numberOfBins,
  //          vtkm::Range& rangeOfValues,
  //          FieldType& binDelta,
  //          vtkm::cont::ArrayHandle<vtkm::Id>& binArray)
  // {
  //   const vtkm::Vec<FieldType, 2> initValue{ vtkm::cont::ArrayGetValue(0, fieldArray) };

  //   vtkm::Vec<FieldType, 2> result =
  //     vtkm::cont::Algorithm::Reduce(fieldArray, initValue, vtkm::MinAndMax<FieldType>());

  //   this->Run(fieldArray, cutInfo, numberOfBins, result[0], result[1], binDelta, binArray);

  //   //update the users data
  //   rangeOfValues = vtkm::Range(result[0], result[1]);
  // }

  template <typename FieldType, typename Storage>
  void PrintArray(vtkm::cont::ArrayHandle<FieldType, Storage> arr)
  {
    int i ;
    for (i = 0; i < arr.GetNumberOfValues(); i++)
    {
      std::cout <<  arr.GetPortalConstControl().Get(i) <<' ';//<<  arr.GetPortalConstControl().Get(i)[1] << " ";
    }
    std::cout << '\n' ;
  }
  // Execute the histogram binning filter given data and number of bins, min,
  // max values.
  // Returns:
  // number of values in each bin
  template <typename FieldType, typename Storage  >
  void Run(vtkm::cont::ArrayHandle<FieldType, Storage> fieldArray,
           vtkm::Vec<vtkm::Id, 3> cutInfo,
           vtkm::Id numberOfBins,
           vtkm::cont::ArrayHandle<vtkm::Float32>& blockMin,
           vtkm::cont::ArrayHandle<vtkm::Float32>& binDelta,
           vtkm::cont::ArrayHandle<vtkm::Id>& binArray)
  {
    const vtkm::Id numberOfValues = fieldArray.GetNumberOfValues();
    vtkm::cont::ArrayHandleCounting<vtkm::Id> eleIndex(0, 1, numberOfValues);
    //std::cout << eleIndex.GetNumberOfValues() << '\n';
    //PrintArray(eleIndex);
    //const auto fieldDelta = compute_delta(fieldMinValue, fieldMaxValue, numberOfBins);

    // Worklet fills in the bin belonging to each value
    vtkm::cont::ArrayHandle<vtkm::Id> blockIndex;
    blockIndex.Allocate(numberOfValues);
    vtkm::cont::ArrayHandle<vtkm::Id> binIndex;
    binIndex.Allocate(numberOfValues);

    // Worklet to set the block number for each data value
    SetBlock blockWorklet(cutInfo);
    vtkm::worklet::DispatcherMapField<SetBlock> setBlockDispatcher(
      blockWorklet);
    setBlockDispatcher.Invoke(eleIndex, blockIndex);

    vtkm::cont::ArrayHandleZip<vtkm::cont::ArrayHandle<vtkm::Id> , vtkm::cont::ArrayHandle<FieldType> > zips = vtkm::cont::make_ArrayHandleZip( blockIndex, fieldArray);
    // Sort the resulting bin array for counting
    vtkm::cont::Algorithm::Sort(zips);
    vtkm::worklet::Zip2 unzip;
    unzip.Run(zips,blockIndex,fieldArray);

    vtkm::cont::ArrayHandle<vtkm::Id> uniqueKeys;
    vtkm::cont::ArrayHandle<FieldType> blockMax;
    vtkm::cont::Algorithm::ReduceByKey(blockIndex, fieldArray, uniqueKeys, blockMax, vtkm::Maximum());
    vtkm::cont::Algorithm::ReduceByKey(blockIndex, fieldArray, uniqueKeys, blockMin, vtkm::Minimum());
    //PrintArray(BlockMax);
    std::cout << "uniquekeys:"<<uniqueKeys.GetNumberOfValues() << '\n';
    //std::cout <<blockMin.GetPortalConstControl().Get(10000) << '\n';
    
    ComputeDelta<FieldType> deltaWorklet(numberOfBins);
    vtkm::worklet::DispatcherMapField<ComputeDelta<FieldType>> computeDeltaDispatcher(deltaWorklet);
    computeDeltaDispatcher.Invoke(blockMin, blockMax, binDelta);
    
    //std::cout <<binDelta.GetPortalConstControl().Get(7999) << '\n';
    // std::cout << "blockIndex:"<<blockIndex.GetNumberOfValues() << '\n';
    // std::cout << "fieldArray:"<<fieldArray.GetNumberOfValues() << '\n';
    // std::cout << "blockIndex:"<<blockIndex.GetNumberOfValues() << '\n';
    // vtkm :: cont :: Algorithm :: Sort ( blockIndex );
    // std::cout <<"Min:" << blockIndex.GetPortalConstControl().Get(0) << '\n';
    // vtkm :: cont :: Algorithm :: Sort ( blockIndex , vtkm :: SortGreater ());
    // std::cout <<"Max:" << blockIndex.GetPortalConstControl().Get(0) << '\n';
    // vtkm::cont::ArrayHandle<FieldType> temp;
    // vtkm::cont::ArrayCopy(blockIndex , temp );
    // std::cout <<"copy:" << temp.GetPortalConstControl().Get(0) << '\n';
    


    
    // Worklet to set the bin number for each data value
    SetHistogramBin<FieldType> binWorklet(numberOfBins,blockMin,binDelta);
    vtkm::worklet::DispatcherMapField<SetHistogramBin<FieldType>> setHistogramBinDispatcher(binWorklet);
    setHistogramBinDispatcher.Invoke(blockIndex, fieldArray, binIndex);
    
    
    //PrintArray(binIndex);
    vtkm::cont::ArrayHandleZip<vtkm::cont::ArrayHandle<vtkm::Id> , vtkm::cont::ArrayHandle<vtkm::Id> > zips2 = vtkm::cont::make_ArrayHandleZip( blockIndex, binIndex);
    // Sort the resulting bin array for counting
    vtkm::cont::Algorithm::Sort(zips2);
    unzip.Run(zips2,blockIndex,binIndex);
    //PrintArray(binIndex);

    BlockShift blockShiftWorklet(numberOfBins);
    vtkm::worklet::DispatcherMapField<BlockShift> BlockShiftkDispatcher(
      blockShiftWorklet);
    BlockShiftkDispatcher.Invoke(blockIndex, binIndex , binIndex);
    
    // Get the upper bound of each bin number
    vtkm::cont::ArrayHandle<vtkm::Id> totalCount;
    vtkm::cont::ArrayHandleCounting<vtkm::Id> binCounter(0, 1, numberOfBins*cutInfo[0]*cutInfo[1]*cutInfo[2]);
    vtkm::cont::Algorithm::UpperBounds(binIndex, binCounter, totalCount);

    // Difference between adjacent items is the bin count
    vtkm::worklet::DispatcherMapField<AdjacentDifference> dispatcher;
    dispatcher.Invoke(binCounter, totalCount, binArray);
    
    //update the users data
    //binDelta = fieldDelta;
    
    std::cout << "FK" << std::endl;

  }
};
}
} // namespace vtkm::worklet

#endif // vtk_m_worklet_FieldHistogram_h
