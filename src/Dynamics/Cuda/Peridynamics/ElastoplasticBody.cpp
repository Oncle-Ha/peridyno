#include "ElastoplasticBody.h"

#include "Topology/TriangleSet.h"
#include "Topology/PointSet.h"
#include "Collision/NeighborPointQuery.h"

#include "Mapping/PointSetToPointSet.h"

#include "Module/Peridynamics.h"
#include "Module/ElastoplasticityModule.h"

#include "Auxiliary/DataSource.h"

#include "ParticleSystem/Module/PositionBasedFluidModel.h"
#include "ParticleSystem/Module/ParticleIntegrator.h"
#include "ParticleSystem/Module/DensityPBD.h"
#include "ParticleSystem/Module/ImplicitViscosity.h"

#include "SharedFunc.h"

namespace dyno
{
	IMPLEMENT_TCLASS(ElastoplasticBody, TDataType)

	template<typename TDataType>
	ElastoplasticBody<TDataType>::ElastoplasticBody()
		: ParticleSystem<TDataType>()
	{
		auto horizon = std::make_shared<FloatingNumber<TDataType>>();
		horizon->varValue()->setValue(Real(0.0085));

		auto m_integrator = std::make_shared<ParticleIntegrator<TDataType>>();
		this->stateTimeStep()->connect(m_integrator->inTimeStep());
		this->statePosition()->connect(m_integrator->inPosition());
		this->stateVelocity()->connect(m_integrator->inVelocity());
		this->stateForce()->connect(m_integrator->inForceDensity());
		this->animationPipeline()->pushModule(m_integrator);
		
		auto m_nbrQuery = std::make_shared<NeighborPointQuery<TDataType>>();
		horizon->outFloating()->connect(m_nbrQuery->inRadius());
		this->statePosition()->connect(m_nbrQuery->inPosition());
		this->animationPipeline()->pushModule(m_nbrQuery);

		auto m_plasticity = std::make_shared<ElastoplasticityModule<TDataType>>();
		horizon->outFloating()->connect(m_plasticity->inHorizon());
		this->stateTimeStep()->connect(m_plasticity->inTimeStep());
		this->statePosition()->connect(m_plasticity->inPosition());
		this->stateVelocity()->connect(m_plasticity->inVelocity());
		this->stateRestShape()->connect(m_plasticity->inRestShape());
		m_nbrQuery->outNeighborIds()->connect(m_plasticity->inNeighborIds());
		this->animationPipeline()->pushModule(m_plasticity);

		auto m_visModule = std::make_shared<ImplicitViscosity<TDataType>>();
		m_visModule->varViscosity()->setValue(Real(1));
		this->stateTimeStep()->connect(m_visModule->inTimeStep());
		horizon->outFloating()->connect(m_visModule->inSmoothingLength());
		this->statePosition()->connect(m_visModule->inPosition());
		this->stateVelocity()->connect(m_visModule->inVelocity());
		m_nbrQuery->outNeighborIds()->connect(m_visModule->inNeighborIds());
		this->animationPipeline()->pushModule(m_visModule);
	}

	template<typename TDataType>
	ElastoplasticBody<TDataType>::~ElastoplasticBody()
	{
		
	}

	template<typename TDataType>
	void ElastoplasticBody<TDataType>::resetStates()
	{
		ParticleSystem<TDataType>::resetStates();

		auto nbrQuery = std::make_shared<NeighborPointQuery<TDataType>>();
		m_horizon.connect(nbrQuery->inRadius());
		this->statePosition()->connect(nbrQuery->inPosition());
		nbrQuery->update();

		if (!this->statePosition()->isEmpty())
		{
			this->stateRestShape()->allocate();
			auto nbrPtr = this->stateRestShape()->getDataPtr();
			nbrPtr->resize(nbrQuery->outNeighborIds()->getData());

			constructRestShape(*nbrPtr, nbrQuery->outNeighborIds()->getData(), this->statePosition()->getData());
		}
	}

	template<typename TDataType>
	void ElastoplasticBody<TDataType>::updateTopology()
	{
		auto ptSet = this->statePointSet()->getDataPtr();
		auto& pts = ptSet->getPoints();
		pts.assign(this->statePosition()->getData());

		auto tMappings = this->getTopologyMappingList();
		for (auto iter = tMappings.begin(); iter != tMappings.end(); iter++)
		{
			(*iter)->apply();
		}
	}

	DEFINE_CLASS(ElastoplasticBody);
}