#include "Cloth.h"
#include "Topology/TriangleSet.h"
#include "Topology/PointSet.h"
#include "Mapping/PointSetToPointSet.h"

#include "ParticleSystem/Module/ParticleIntegrator.h"

#include "Collision/NeighborPointQuery.h"

#include "Module/LinearElasticitySolver.h"
#include "Module/Peridynamics.h"
#include "Module/FixedPoints.h"

#include "Auxiliary/DataSource.h"

#include "SharedFunc.h"
#include "TriangularSystem.h"

namespace dyno
{
	IMPLEMENT_TCLASS(Cloth, TDataType)

	template<typename TDataType>
	Cloth<TDataType>::Cloth()
		: TriangularSystem<TDataType>()
	{
		auto horizon = std::make_shared<FloatingNumber<TDataType>>();
		horizon->varValue()->setValue(Real(0.01));
		this->animationPipeline()->pushModule(horizon);

		auto integrator = std::make_shared<ParticleIntegrator<TDataType>>();
		this->stateTimeStep()->connect(integrator->inTimeStep());
		this->statePosition()->connect(integrator->inPosition());
		this->stateVelocity()->connect(integrator->inVelocity());
		this->stateForce()->connect(integrator->inForceDensity());

		this->animationPipeline()->pushModule(integrator);

		auto elasticity = std::make_shared<LinearElasticitySolver<TDataType>>();
		horizon->outFloating()->connect(elasticity->inHorizon());
		this->stateTimeStep()->connect(elasticity->inTimeStep());
		this->statePosition()->connect(elasticity->inPosition());
		this->stateVelocity()->connect(elasticity->inVelocity());
		this->stateRestShape()->connect(elasticity->inRestShape());
		this->animationPipeline()->pushModule(elasticity);
	}

	template<typename TDataType>
	Cloth<TDataType>::~Cloth()
	{
		
	}

	template<typename TDataType>
	void Cloth<TDataType>::resetStates()
	{
		auto input = this->inTriangleSet()->getDataPtr();

		auto ts = this->stateTriangleSet()->getDataPtr();

		ts->copyFrom(*input);

		//auto& nbr = nbrQuery->outNeighborIds()->getData();
		auto& nbr = ts->getPointNeighbors();
		auto& pts = ts->getPoints();

		this->statePosition()->resize(pts.size());
		this->stateVelocity()->resize(pts.size());
		this->stateForce()->resize(pts.size());

		this->statePosition()->assign(pts);
		this->stateVelocity()->reset();

		if (this->stateRestShape()->isEmpty())
		{
			this->stateRestShape()->allocate();
		}
		
		auto nbrPtr = this->stateRestShape()->getDataPtr();

		this->stateOldPosition()->assign(this->statePosition()->getData());

		constructRestShapeWithSelf(*nbrPtr, nbr, this->statePosition()->getData());
	}

	template<typename TDataType>
	void Cloth<TDataType>::preUpdateStates()
	{
		auto& posOld = this->stateOldPosition()->getData();
		auto& posNew = this->statePosition()->getData();

		posOld.assign(posNew);
	}

	DEFINE_CLASS(Cloth);
}