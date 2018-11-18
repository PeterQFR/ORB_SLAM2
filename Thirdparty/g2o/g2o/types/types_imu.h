// g2o - General Graph Optimization
// Copyright (C) 2011 H. Strasdat
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// * Redistributions of source code must retain the above copyright notice,
//   this list of conditions and the following disclaimer.
// * Redistributions in binary form must reproduce the above copyright
//   notice, this list of conditions and the following disclaimer in the
//   documentation and/or other materials provided with the distribution.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
// IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
// TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
// PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
// TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

// Modified by Raúl Mur Artal (2014)
// Added EdgeSE3ProjectXYZ (project using focal_length in x,y directions)
// Modified by Raúl Mur Artal (2016)
// Added EdgeStereoSE3ProjectXYZ (project using focal_length in x,y directions)
// Added EdgeSE3ProjectXYZOnlyPose (unary edge to optimize only the camera pose)
// Added EdgeStereoSE3ProjectXYZOnlyPose (unary edge to optimize only the camera pose)

#ifndef G2O_SIX_DOF_TYPES_EXPMAP
#define G2O_SIX_DOF_TYPES_EXPMAP

#include "../core/base_vertex.h"
#include "../core/base_binary_edge.h"
#include "../core/base_unary_edge.h"
#include "se3_ops.h"
#include "se3quat.h"
#include <si3quat.h>
#include "types_sba.h"
#include <Eigen/Geometry>

namespace g2o {
namespace types_imu {
void init();
}

using namespace Eigen;

typedef Matrix<double, 6, 6> Matrix6d;
typedef Vector<doulbe, 1, 15> Vector15d;

/**
 * \brief SE3 Vertex parameterized internally with a transformation matrix
 and externally with its exponential map
 */
class  VertexSI3 : public BaseVertex<16, SI3Quat>{
public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  VertexSI3();

  bool read(std::istream& is);

  bool write(std::ostream& os) const;

  virtual void setToOriginImpl() {
    _estimate._t = SE3Quat();
    _estimate._ba = Vector3d::Zero();
    _estimate._bw = Vector3d::Zero();
    _estimate._v = Vector3d::Zero();
  }

  //Provide an update of the Vertex by applying a small movement in various
  //Direction of type SI3. This Update is the minimal representation
  //ie quaternion axis and vector3 for velocity and position.
  //Update for bias is a 3d vector for each accleration and gyro.
  virtual void oplusImpl(const double* update_)  {

	Eigen::Map<const Vector15d> update(update_);
	SE3Quat se3Q = SE3Quat::exp(update.block(0,0,1,6))*estimate()._t;

	Vector3d v = estimate()._v + update.block(0, 6, 1,3);
	Vector3d ba = estimate()._ba + update.block(0, 9, 1, 3);
	Vector3d bw = estimate()._bw + update.block(0, 12, 1, 3);

	SI3Quat si3(se3Q, v, ba, bw);

	setEstimate(si3);

  }

};

class VertexVelocity : public BaseVertex<3, Vector3d>
{
public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW

	VertexVelocity():BaseVertex<3, Vector3d>(){};

	bool read(std::istream& is);
	bool write(Std::ostream& os) const;

	virtual void setToOriginalImple(){
		_estimate = Vector3d::Zero();
	}

	virtual void oplusImpl (const double* update_){
		Eigen::Map<const Vector3d> update(update_);

		setEstimate(estimate() + update);
	}
};

class VertexBias : public BaseVertex<6, Vector6d>
{
public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW

	VertexBias():BaseVertex<6, Vector6d>(){};

	bool read(std::istream& is);
	bool write(Std::ostream& os) const;

	virtual void setToOriginalImple(){
		_estimate = Vector6d::Zero();
	}

	virtual void oplusImpl (const double* update_){
		Eigen::Map<const Vector6d> update(update_);

		setEstimate(estimate() + update);
	}

	Vector3d getGyroBias(){
		return _estimate.block(0,3,1,3);
	}

	Vector3d getAccelBias(){
		return _estimate.block(0, 0, 1, 3);
	}

};

/**
 * The Measurement that is used is the preintegral from imu data
 * containing Pose, Velocity and Bias.
 */
class EdgeImuUpdate : public BaseMultiEdge<6, SI3Quat>
{
public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW

	EdgeImuUpdate(): BaseMultiEdge<15, SI3Quat>()
	{
		//Verticies are a specific size. 6 verticies are added.
		/**
		 * They are:
		 * Pi,
		 * vi,
		 * bi,
		 * Pj,
		 * vj,
		 * bj
		 */
		_vertices.resize(6);
	};
	/**
	 * \brief _measurement is the preintegral.
	 */
	void computeError()  {
		//here we determine the er

	    //The error is given by 6 in R.Mur-artal 2016. visual inertial monocular
	    //SLAM with Map reuse.
		SE3Quat xi = static_cast<const VertexSE3Expmap*>(
				_vertices[0])->estimate();
		Vector3d vi = static_cast<const VertexVelocity*>(_vertices[1])->estimate();
		Vector3d bgi = static_cast<const VertexBias*>(_vertices[2])->getGyroBias();
		Vector3d bai = static_cast<const VertexBias*>(_vertices[2])->getAccelBias();

		SE3Quat xj = static_cast<const VertexSE3Expmap*>(
						_vertices[3])->estimate();
		Vector3d vj = static_cast<const VertexVelocity*>(_vertices[4])->estimate();
		Vector3d bgj = static_cast<const VertexBias*>(_vertices[5])->getGyroBias();
		Vector3d baj = static_cast<const VertexBias*>(_vertices[5])->getAccelBias();
		Vector3d g(0.0,0.0,-9.8);

		//preintegral
		SI3Quat preIntegral = this->_measurement;
		Matrix3d Ribw = xi._r.conjugate().toRotationMatrix();
		Matrix3d Rjwb = xj._r.toRotationMatrix();
		Matrix3d dRot = preIntegral._t._r.toRotationMatrix();
		Vector3d eR = log((dRot*exp(_gJdR*bgj)).transpose()*
				Ribw*Rjwb);
		Vector3d eV = Ribw *(vj - vi - g*dt) - (preIntegral._v +
				_gJdv*bgj + _aJdv*baj);
		Vector3d eP = Ribw(xj._t - xi._t - _vi*dt - 0.5*g*dt*dt) -
				(preIntegral._t._t + _gJdp*bgj + _aJdp*baj);
		Vector6d eB = static_cast<const VertexBias*>(_vertices[2])->getEstimateData() -
				static_cast<const VertexBias*>(_vertices[5])->getEstimateData();

		SI3Quat residual(Quaterniond(skew(eR)), eP, eV, eB.block(0,0,1,3),
				eB.block(0,3,1,3));

	}

	/**
	 * \brief preIntegrate takes in accleration and gyro data and evaluates
	 * a preintegral which is of type SI3Quat.
	 * \param accel is a vector of accel measuremnts in m/S^2)
	 * \param gyro is a vector of gyro measuremetns in rad/s
	 * \param timeoffsets of the above measurements from the anchor time in seconds
	 */
	SI3Quat calculatePreIntegral(const std::vector<Vector3d>& accel,
			const std::vector<Vector3d>& gyro,
			const std::vector<double>& timeoffsets)
	{
		const VertexSE3Expmap* xi = static_cast<const VertexSE3Expmap*>(
				_vertices[0]);
		const VertexBias* bi = static_cast<const VertexBias*>(_vertices[2]);
		const VertexVelocity* vi = static_cast<const VertexVelocity*> (_vertices[1]);

		SI3Quat anchor(xi->estimate()._r, xi->estimate()._t, vi->estimate(),
				bi->estimate().getAccelBias(), bi->estimate().getGyroBias());

		for (int k =0; k< accel.size(); k++)
		{
			double dt;
			if(k==0)
				dt = timeoffsets[k];
			else
				dt = timeoffsets[k]-timeoffsets[k-1];
			anchor.preIntegrate(accel[k], gyro[k], dt);

			//Here is a good place to now calculate the Jacobians for this
			//step.
			Matrix3d Rik = anchor._t._r.toRotationMatrix();
			_gJdR += -(Rik*bi->getGyroBias()
					*rightJacobian(gyro[k]*0.01)*dt);

			Matrix3d accelSkew = Rik*skew(accel[k] -bi->getAccelBias()) *_gJdR*dt;

			_aJdv += -(Rik*dt);
			_gJdv += -(accelSkew);
			_aJdp += _aJdv*dt - 0.5*Rik*dt*dt;
			_gJdp += _gJdv*dt - 0.5*accelSkew*dt;

		}

		this->setMeasurement(anchor);
	}
protected:
	Vector3d& _aJdp;
	Vector3d& _gJdp;
	Vector3d& _gJdv;
	Vector3d& _aJdv;
	Vector3d& _gJdR;
public:

	/**
	 * \brief as this is a multi constraint. The ordering of the vertices are
	 * essential to ensure all information is correctly allocated.
	 */
	/**
	 * \brief sets the prior pose P. This is from a previous iteration or from
	 * a keyframe if a map update has just occurred.
	 */
	inline void setPoseIVertex(const VertexSE3Expmap* vertex){
		setVertex(0, dynamic_cast<g2o::OptimizableGraph::Vertex*>(vertex));
	}
	inline void setVelocityIVertex(const VertexVelocity* vertex)
	{
		setVertex(1, dynamic_cast<g2o::OptimizableGraph::Vertex*>(vertex));
	}

	inline void setBiasIVertex(const VertexBias* vertex)
	{
		setVertex(2, dynamic_cast<g2o::OptimizableGraph::Vertex*>(vertex));
	}


	inline void setPoseJVertex( const VertexSE3Expmap* vertex){
		setVertex(3, dynamic_cast<g2o::OptimizableGraph::Vertex*>(vertex));
	}
	inline void setVelocityJVertex(const VertexVelocity* vertex)
	{
		setVertex(4, dynamic_cast<g2o::OptimizableGraph::Vertex*>(vertex));
	}

	inline void setBiasJVertex(const VertexBias* vertex)
	{
		setVertex(5, dynamic_cast<g2o::OptimizableGraph::Vertex*>(vertex));
	}

protected:



	inline Vector3d calculateRotResidual(const Matrix3d& preRot, const Vector3d& wdt,
			const Vector3d& bgyroj, const Matrix3d& roti, const Matrix3d* rotj)
	{
		return log((preRot*exp())
	}

};

class EdgeBiasUpdate : public BaseBinaryEdge<6, Vector6d, VertexBias, VertexBias>
{
public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW

	EdgeBiasUpdate():BaseBinaryEdge<6, Vector6d, VertexBias, VertexBias>(){};
	bool read(std::istream& is);
	bool write(std::ostream& os) const;

	void computeError(){
		const VertexBias* bj = static_cast<const VertexBias*>(_vertices[1]);
	    const VertexBias* bi = static_cast<const VertexBias*>(_vertices[0]);
	    _error = bj - bi ;
	};


} // end namespace

#endif
