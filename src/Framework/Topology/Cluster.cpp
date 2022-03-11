#pragma once
#include "Cluster.h"

namespace dyno
{
    // IMPLEMENT_TCLASS(Cluster, TDataType)
    template<typename TDataType>
    Cluster<TDataType>::Cluster()
    {
    }
    
	template<typename TDataType>
    Cluster<TDataType>::Cluster(
        const int* indices, int indicesCount, 
        const double* weights, int weightsCount, 
        double* Mt, double* Mtl)
    {
        std::vector<int> tmp_indices(indices, indices + indicesCount);
        m_indices.resize(indicesCount);
		m_indices.assign(tmp_indices);

        std::vector<float> tmp_weights(weights, weights + weightsCount);
        m_weights.resize(weightsCount);
        m_weights.assign(tmp_weights);

        m_transform = Mat4f(Mt[0], Mt[1], Mt[2], Mt[3],
                            Mt[4], Mt[5], Mt[6], Mt[7],
                            Mt[8], Mt[9], Mt[10], Mt[11],
                            Mt[12], Mt[13], Mt[14], Mt[15]);

        m_transformLink = Mat4f(Mtl[0], Mtl[1], Mtl[2], Mtl[3],
                            Mtl[4], Mtl[5], Mtl[6], Mtl[7],
                            Mtl[8], Mtl[9], Mtl[10], Mtl[11],
                            Mtl[12], Mtl[13], Mtl[14], Mtl[15]);                        
    }

    template<typename TDataType>
    Cluster<TDataType>::~Cluster()
    {

    }
#ifdef PRECISION_FLOAT
	template class Cluster<DataType3f>;
#else
	template class Cluster<DataType3d>;
#endif    
}