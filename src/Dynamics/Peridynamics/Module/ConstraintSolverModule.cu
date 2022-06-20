#include "ConstraintSolverModule.h"
#include "Node.h"
#include "Matrix/MatrixFunc.h"
#include "ParticleSystem/Kernel.h"

namespace dyno
{
	IMPLEMENT_TCLASS(ConstraintSolverModule, TDataType)
	
	template<typename TDataType>
	ConstraintSolverModule<TDataType>::ConstraintSolverModule()
		: ConstraintModule()
	{
		this->inHorizon()->setValue(0.0125);
		this->inNeighborIds()->tagOptional(true);
	}


	template<typename TDataType>
	ConstraintSolverModule<TDataType>::~ConstraintSolverModule()
	{
        // clear()
	}


	template<typename TDataType>
	void ConstraintSolverModule<TDataType>::constrain()
	{
        // solve()
	}


	template<typename TDataType>
	void ConstraintSolverModule<TDataType>::preprocess()
	{
        //resize 
		//this->computeMaterialStiffness();
	}

	DEFINE_CLASS(ConstraintSolverModule);
}