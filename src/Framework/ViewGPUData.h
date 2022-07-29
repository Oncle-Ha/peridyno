    #pragma once
#include "Vector.h"
#include "Module.h"
#include <Eigen/Core>
#include "Module/TopologyModule.h"

#define LIMIT_SIZE 10000
#define LINE_SIZE 100
#define ROW_SIZE 100
namespace dyno
{
	/*!
	*	\class	ViewGPUData
	*	\brief	View GPU Data with Image Watch
	*/    
	template<typename TDataType>
	class ViewGPUData
	{
	public:

		typedef typename TDataType::Real Real;
		typedef typename TDataType::Coord Coord;
		typedef typename TDataType::Matrix Matrix;
		typedef typename TopologyModule::Pair2 Pair2;
        typedef typename TopologyModule::Pair3 Pair3;		
        typedef typename Eigen::Matrix<Real, LINE_SIZE, ROW_SIZE> EMat1000x10;
        
		typedef VectorND<Real, 3> Pair3f;

		ViewGPUData();
		ViewGPUData(int size);
		~ViewGPUData();

        DArray<Real>& getGPU() {return m_GPU;}
        CArray<Real>& getCPU() {return m_CPU;}
        
        void resize(int size);
		
		bool viewIm(DArray<Real> &p_GPU);
		bool viewIm(DArray<Coord> &p_GPU, int index);
		bool viewIm(DArray<Matrix> &p_GPU, int row, int line);
		bool viewIm(DArray<Pair2> &p_GPU, int index);
		bool viewIm(DArray<Pair3> &p_GPU, int index);
		bool viewIm(DArray<Pair3f> &p_GPU, int index);

        bool view();

	private:
		DArray<Real> m_GPU;
        CArray<Real> m_CPU;
        int m_size;
        EMat1000x10 m_view;
	};
}