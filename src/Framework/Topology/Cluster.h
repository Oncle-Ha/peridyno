#pragma once
#include "Module.h"

namespace dyno
{
	/*!
	*	\class	Cluster
	*	\brief	Each cluster acts on a subset of the geometry's control points, with different weights. 
	*/    
    template<typename TDataType>
	class Cluster
	{ 
        // DECLARE_CLASS_1(Cluster, TDataType)
	public:
		typedef typename TDataType::Real Real;
		typedef typename TDataType::Coord Coord;
		typedef typename TDataType::Matrix Matrix;

        Cluster();
        Cluster(const int*, int, const double*, int, double*, double*);
        ~Cluster();
        

        DArray<int> m_indices;
        DArray<Real> m_weights;
        // Skeleton -> Mesh : M_t * M_tl^-1
        Mat4f m_transform; // M_t ?
        Mat4f m_transformLink; // M_tl

        int m_jointIndex; // m_jointMap[index]
    };
}