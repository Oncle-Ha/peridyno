#pragma once
#include "Cluster.h"

namespace dyno
{
    template<typename TDataType>
    Cluster<TDataType>::Cluster()
    {
    }
    

#ifdef PRECISION_FLOAT
	template class Cluster<DataType3f>;
#else
	template class Cluster<DataType3d>;
#endif    
}