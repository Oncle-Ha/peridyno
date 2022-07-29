#pragma once
#include "ViewGPUData.h"

namespace dyno
{
    template<typename TDataType>
    ViewGPUData<TDataType>::ViewGPUData()
    {
        m_view = EMat1000x10::Zero();
        m_size = 0;
    }
    
    template<typename TDataType>
    ViewGPUData<TDataType>::ViewGPUData(int size)
    {
        m_view = EMat1000x10::Zero();
        this->resize(size);
    }

    template<typename TDataType>
    ViewGPUData<TDataType>::~ViewGPUData()
    {
        m_GPU.clear();
        m_CPU.clear();
    }

    template<typename TDataType>
    void ViewGPUData<TDataType>::resize(int size)
    {
        if (size > LIMIT_SIZE) size = LIMIT_SIZE;
        m_GPU.resize(size);
        m_size = size;
    }

    template<typename TDataType>
    bool ViewGPUData<TDataType>::view()
    {
        if (m_size != m_GPU.size()) return false;

        m_CPU.assign(m_GPU);
        for (int i = 0; i < m_size; ++i)
        {
            m_view(i % LINE_SIZE, i / LINE_SIZE) = m_CPU[i];
        }

        return true;
    }


    template<typename TDataType>
    bool ViewGPUData<TDataType>::viewIm(DArray<Real> &p_GPU)
    {
        m_GPU.assign(p_GPU, m_size, 0, 0);
        this->view();
        return true;
    }    

    template <typename Coord, typename Real>
	__global__ void VGPUCoord(
        DArray<Coord> p,
        DArray<Real> m,
        int index,
        int size)
    {
        int pId = threadIdx.x + (blockIdx.x * blockDim.x);
		if(pId >= size) return;    
        m[pId] = p[pId][index];
    }

    template<typename TDataType>
    bool ViewGPUData<TDataType>::viewIm(DArray<Coord> &p_GPU, int index)
    {
        cuExecute(m_size,
            VGPUCoord,
            p_GPU,
            m_GPU,
            index,
            m_size);
        cuSynchronize();
        this->view();
        return true;
    }   


    template <typename Matrix, typename Real>
	__global__ void VGPUMatrix(
        DArray<Matrix> p,
        DArray<Real> m,
        int row,
        int line,
        int size)
    {
        int pId = threadIdx.x + (blockIdx.x * blockDim.x);
		if(pId >= size) return;    
        m[pId] = p[pId](row, line);
    }

    template<typename TDataType>
    bool ViewGPUData<TDataType>::viewIm(DArray<Matrix> &p_GPU, int row, int line)
    {
        cuExecute(m_size,
            VGPUMatrix,
            p_GPU,
            m_GPU,
            row,
            line,
            m_size);
        cuSynchronize();
        this->view();        
        return true;
    }   

    template <typename Pair2, typename Real>
	__global__ void VGPUPair2(
        DArray<Pair2> p,
        DArray<Real> m,
        int index,
        int size)
    {
        int pId = threadIdx.x + (blockIdx.x * blockDim.x);
		if(pId >= size) return;    
        m[pId] = Real(p[pId][index]);
    }

    template<typename TDataType>
    bool ViewGPUData<TDataType>::viewIm(DArray<Pair2> &p_GPU, int index)
    {
        cuExecute(m_size,
            VGPUPair2,
            p_GPU,
            m_GPU,
            index,
            m_size);
        cuSynchronize();
        this->view();         
        return true;
    }   

    template <typename Pair3, typename Real>
	__global__ void VGPUPair3(
        DArray<Pair3> p,
        DArray<Real> m,
        int index,
        int size)
    {
        int pId = threadIdx.x + (blockIdx.x * blockDim.x);
		if(pId >= size) return;    
        m[pId] = Real(p[pId][index]);
    }

    template<typename TDataType>
    bool ViewGPUData<TDataType>::viewIm(DArray<Pair3> &p_GPU, int index)
    {
        cuExecute(m_size,
            VGPUPair3,
            p_GPU,
            m_GPU,
            index,
            m_size);
        cuSynchronize();
        this->view();          
        return true;
    }       

    template <typename Pair3f, typename Real>
	__global__ void VGPUPair3f(
        DArray<Pair3f> p,
        DArray<Real> m,
        int index,
        int size)
    {
        int pId = threadIdx.x + (blockIdx.x * blockDim.x);
		if(pId >= size) return;    
        m[pId] = Real(p[pId][index]);
    }

    template<typename TDataType>
    bool ViewGPUData<TDataType>::viewIm(DArray<Pair3f> &p_GPU, int index)
    {
        cuExecute(m_size,
            VGPUPair3f,
            p_GPU,
            m_GPU,
            index,
            m_size);
        cuSynchronize();
        this->view();
        return true;
    }            

#ifdef PRECISION_FLOAT
	template class ViewGPUData<DataType3f>;

#else
    template class ViewGPUData<DataType3d>;
#endif    
}