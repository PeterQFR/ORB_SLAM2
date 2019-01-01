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

#ifndef G2O_TYPES_IMU
#define G2O_TYPES_IMU

#include "../core/base_vertex.h"
#include "../core/base_binary_edge.h"
#include "../core/base_unary_edge.h"
#include "../core/base_multi_edge.h"
#include "se3_ops.h"
#include "se3quat.h"
#include "si3quat.h"
#include "types_sba.h"
#include "Eigen/Geometry"
#include "types_six_dof_expmap.h"


namespace g2o {
namespace types_imu {
void init();
}

using namespace Eigen;

typedef Matrix<double, 6, 6> Matrix6d;
typedef Matrix<double, 1, 15> Vector15d;
typedef Matrix<double, 9, 9> Matrix9d;
typedef Matrix<double, 9, 6> ControlMatrix;


struct StateInfo {
	g2o::SE3Quat pose;
	Vector3d vel;
	Vector3d bg;
	Vector3d ba;
	double stamp;
};

struct InertialMeasurement {
	Vector3d gyro;
	Vector3d acceleration;
	double stamp;
};


class VertexVelocity : public BaseVertex<3, Vector3d>
{
public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW

	VertexVelocity(): BaseVertex<3, Vector3d>(){};

	bool read(std::istream& is);
	bool write(std::ostream& os) const;

	virtual void setToOriginImpl(){
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
	bool write(std::ostream& os) const;


	virtual void setToOriginImpl(){
		_estimate = Vector6d::Zero();
	};

	virtual void oplusImpl(const double* update_){
		Eigen::Map<const Vector6d> update(update_);
		setEstimate(estimate() + update);
	};
public:
	Vector3d getGyroBias( ) const {
		return _estimate.tail(3);
	}

	Vector3d getAccelBias() const {
		return _estimate.head(3);
	}

	void setGyroBias(const Vector3d& bias){
		_estimate.tail(3) = bias.head(3);
	}
	void setAccelBias(const Vector3d& bias)
	{
		_estimate.head(3)=bias.head(3);
	}


};

/**
 * The Measurement that is used is the preintegral from imu data
 * containing Pose, Velocity and Bias.
 */
class EdgeImuUpdate : public BaseMultiEdge<15, SI3Quat>
{
public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW

	EdgeImuUpdate() :
	BaseMultiEdge<15, SI3Quat>(),
	_aJdp(Matrix3d::Zero()),
	_gJdp(Matrix3d::Zero()),
	_gJdv(Matrix3d::Zero()),
	_aJdv(Matrix3d::Zero()),
	_gJdR(Matrix3d::Zero()),
	_dt(0.0),
	_na(0.01),
	_ng(0.001),
	_nba(0.01),
	_nbg(0.01),
	_zetaNoise(Matrix6d::Zero())
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
		resize(6);
		setNoises(_ng, _na, _nba, _nbg);
	};

	void setNoises(const double& ng, const double& na,
			const double& nba, const double& nbg)
	{
		_zetaNoise.block(0,0,3,3) =ng*Matrix3d::Identity();
		_zetaNoise.block(3,3,3,3) = na*Matrix3d::Identity();


		this->_information.block(9,9,3,3) = (1.0/nba)*Matrix3d::Identity();
		this->_information.block(12,12,3,3) = (1.0/nbg)*Matrix3d::Identity();

	}

	bool read(std::istream& is);
	bool write(std::ostream& os) const;
	/**
	 * \brief _measurement is the preintegral.
	 */
	void computeError()  {
		//here we determine the er

	    //The error is given by 6 in R.Mur-artal 2016. visual inertial monocular
	    //SLAM with Map reuse.
		SE3Quat xi = static_cast<const VertexSE3Expmap*>(_vertices[0])->estimate();
		Vector3d vi = static_cast<const VertexVelocity*>(_vertices[1])->estimate();

		SE3Quat xj = static_cast<const VertexSE3Expmap*>(_vertices[3])->estimate();
		Vector3d vj = static_cast<const VertexVelocity*>(_vertices[4])->estimate();
		Vector3d bgj = static_cast<const VertexBias*>(_vertices[5])->getGyroBias();
		Vector3d baj = static_cast<const VertexBias*>(_vertices[5])->getAccelBias();
		Vector3d g(0.0,0.0,-9.8);

		//preintegral
		SI3Quat preIntegral = this->_measurement;
		Matrix3d Ribw = xi.rotation().conjugate().toRotationMatrix();
		Matrix3d Rjwb = xj.rotation().toRotationMatrix();
		Matrix3d dRot = preIntegral.rotation().toRotationMatrix();
		Vector3d eR = log((dRot*exp(_gJdR*bgj)).transpose()*Ribw*Rjwb);
		Vector3d eV = Ribw *(vj - vi - g*_dt) - (preIntegral.velocity() +
				_gJdv*bgj + _aJdv*baj);
		Vector3d eP = Ribw*(xj.translation() - xi.translation() - vi*_dt - 0.5*g*_dt*_dt) -
				(preIntegral.translation() + _gJdp*bgj + _aJdp*baj);
		Vector6d eB = static_cast<const VertexBias*>(_vertices[2])->estimate() -
				static_cast<const VertexBias*>(_vertices[5])->estimate();

		//SI3Quat residual(Quaterniond(skew(eR)), eP, eV, eB.block(0,0,1,3),
		//		eB.block(0,3,1,3));
		_error[0] = eR[0];
		_error[1] = eR[1];
		_error[2] = eR[2];

		_error[3] = eP[0];
		_error[4] = eP[1];
		_error[5] = eP[2];

		_error[6] = eV[0];
		_error[7] = eV[1];
		_error[8] = eV[2];

		_error[9] = eB[0];
		_error[10] = eB[1];
		_error[11] = eB[2];
		_error[12] = eB[3];
		_error[13] = eB[4];
		_error[14] = eB[5];
	}

	/**
	 * \brief preIntegrate takes in accleration and gyro data and evaluates
	 * a preintegral which is of type SI3Quat.
	 * \param accel is a vector of accel measuremnts in m/S^2)
	 * \param gyro is a vector of gyro measuremetns in rad/s
	 * \param dts the times the measurement is active from the previous measure
	 * ment. The first measurement is the time since the anchor time.
	 */
	void calculatePreIntegral(const std::vector<Vector3d>& accel,
			const std::vector<Vector3d>& gyro,
			const std::vector<double>& dts)
	{
		assert(dts.size()>1);

		const VertexSE3Expmap* xi = static_cast<const VertexSE3Expmap*>(
				_vertices[0]);
		const VertexBias* bi = static_cast<const VertexBias*>(_vertices[2]);
		const VertexVelocity* vi = static_cast<const VertexVelocity*> (
				_vertices[1]);

		SI3Quat anchor(xi->estimate().rotation(), xi->estimate().translation(),
				vi->estimate(),
				bi->getAccelBias(),
				bi->getGyroBias());

		Matrix9d covariance = Matrix9d::Zero();
		/**
		 * the times are times that the measurement is from the
		 */
		for (unsigned int k =0; k< accel.size(); k++)
		{
			_dt += dts[k];
			double dt = dts[k];
			anchor.preIntegrate(accel[k], gyro[k], dt);


			Matrix3d Rik = anchor.rotation().toRotationMatrix();

			//Calculate the processNoise:
			Matrix3d integrator = Rik*skew(accel[k] -bi->getAccelBias()) *dt;

			covariance = covarianceMatrix(covariance,
					calculateAMatrix(Rik, integrator, dt),
					calculateBMatrix(Rik, rightJacobian(gyro[k]*0.01), dt));


			//Here is a good place to now calculate the Jacobians for this
			//step.

			_gJdR += -(Rik.transpose()*rightJacobian(gyro[k]*0.01)*dt);

			Matrix3d accelSkew = integrator*_gJdR;

			_aJdv += -(Rik*dt);
			_gJdv += -(accelSkew);
			_aJdp += _aJdv*dt - 0.5*Rik*dt*dt;
			_gJdp += _gJdv*dt - 0.5*accelSkew*dt;

		}

		this->_information.block(0,0,9,9) = covariance.inverse();
		this->setMeasurement(anchor);

	};

	/**
	 * \brief calculate the simple information matrix based on appendix A of
	 * D.Scarramuzza 2016, on manifold preintegration of real time visual inertial
	 * odometry
	 */
	Matrix9d covarianceMatrix(const Matrix9d& zetabar1,
			const Matrix9d& Amatrix,
			const ControlMatrix& Bmatrix)
	{
		return Amatrix*zetabar1*Amatrix.transpose() +
				Bmatrix*_zetaNoise*Bmatrix.transpose();
	};


protected:
	Matrix9d calculateAMatrix(const Matrix3d& dR,
			const Matrix3d integrator,
			const double& dt)
	{
		Matrix9d Amatrix=Matrix9d::Identity();
		//Matrix3d integrator = dR*skew(a-ba)*dt;

		Amatrix.block(0,0,3,3) = dR.transpose();
		Amatrix.block(3,0,3,3) = -integrator;
		Amatrix.block(6,0,3,3) = -0.5*integrator*dt;
		Amatrix.block(6,3,3,3) = dt*Matrix3d::Identity();
		return Amatrix;


	};

	ControlMatrix calculateBMatrix(
			const Matrix3d& dR,
			const Matrix3d& Jr,
			const double& dt)
	{
		ControlMatrix Bmatrix = ControlMatrix::Zero();
		Matrix3d integral = dR*dt;
		Bmatrix.block(0,0,3,3) = Jr*dt;
		Bmatrix.block(3,3,3,3) = integral;
		Bmatrix.block(6,3,3,3) = 0.5*dt*integral;

		return Bmatrix;

	}

protected:
	Matrix3d _aJdp;
	Matrix3d _gJdp;
	Matrix3d _gJdv;
	Matrix3d _aJdv;
	Matrix3d _gJdR;
	double _dt;

	//Acceleration and gyro noises:
	double _na, _ng, _nba, _nbg;
	Matrix6d _zetaNoise;
public:

	/**
	 * \brief as this is a multi constraint. The ordering of the vertices are
	 * essential to ensure all information is correctly allocated.
	 */
	/**
	 * \brief sets the prior pose P. This is from a previous iteration or from
	 * a keyframe if a map update has just occurred.
	 */
	inline void setPoseIVertex(VertexSE3Expmap* vertex){
		setVertex(0, dynamic_cast<g2o::OptimizableGraph::Vertex*>(vertex));
	};
	inline void setVelocityIVertex( VertexVelocity* vertex)
	{
		setVertex(1, dynamic_cast<g2o::OptimizableGraph::Vertex*>(vertex));
	};

	inline void setBiasIVertex(VertexBias* vertex)
	{
		setVertex(2, dynamic_cast<g2o::OptimizableGraph::Vertex*>(vertex));
	};


	inline void setPoseJVertex(VertexSE3Expmap* vertex){
		setVertex(3, dynamic_cast<g2o::OptimizableGraph::Vertex*>(vertex));
	}
	inline void setVelocityJVertex(VertexVelocity* vertex)
	{
		setVertex(4, dynamic_cast<g2o::OptimizableGraph::Vertex*>(vertex));
	};

	inline void setBiasJVertex(VertexBias* vertex)
	{
		setVertex(5, dynamic_cast<g2o::OptimizableGraph::Vertex*>(vertex));
	};

protected:



	/*inline Vector3d calculateRotResidual(const Matrix3d& preRot, const Vector3d& wdt,
			const Vector3d& bgyroj, const Matrix3d& roti, const Matrix3d* rotj)
	{
		return log(preRot*exp());
	};*/

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
	    _error = bj->estimate() - bi->estimate() ;
	};

	void applyInformation(const double& nba, const double& nbg)
	{
		this->_information.block(0,0,3,3) = (1.0/nba)*Matrix3d::Identity();
		this->_information.block(3,3,3,3) = (1.0/nbg)*Matrix3d::Identity();
	}
};

} // end namespace

#endif
