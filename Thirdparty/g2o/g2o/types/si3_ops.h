/*
 * si3_ops.h
 *
 *  Created on: 2Nov.,2018
 *      Author: peter
 */

#ifndef ORB_SLAM2_THIRDPARTY_G2O_G2O_TYPES_SI3_OPS_H_
#define ORB_SLAM2_THIRDPARTY_G2O_G2O_TYPES_SI3_OPS_H_

#include "se3_ops.h"

namespace g2o {
	using namespace Eigen;

	inline Matrix4d omega(const Vector3d& w)
	{
		Matrix3d sk = skew(w);

		Matrix4d omega = Matrix4d::Zero();
		omega.block(0,0,3,3) = -sk;
		omega.block(0,3,3,1) = w;
		omega.block(3,0,1,3) = -w;
		return omega;
	}

	inline Matrix4d zeroOrderIntegrator(const Vector3d& w,
			const double& dt)
	{
		return Matrix4d::Identity() + dt*0.5*omega(w);
	}


	inline Matrix3d exp(const Vector3d& wdt)
	{
		double length = wdt.norm();
		Vector3d unit = wdt/length;
		Matrix3d skewMatrix = skew(unit);
		return Matrix3d::Identity()*cos(length) + sin(length) * skewMatrix +
				(1- cos(length))*skewMatrix*skewMatrix.transpose();
	}

	inline Vector3d log(const Matrix3d& R)
	{
		double phi = acos((R.trace() -1)/2);
		Vector3d unit = Vector3d::Zero();
		unit(0)=R(2,1)-R(1,2);
		unit(1)=R(0,2)-R(2,0);
		unit(2)=R(1,0)-R(0,1);
		unit /= 2*sin(phi);
		return unit*phi;
	}

	inline Matrix3d rightJacobian(const Vector3d& w)
	  {

		  double norm = w.norm();
		  double invnorm = 1.0/norm;
		  Matrix3d sk=skew(w);
		  return Matrix3d::Identity() - (1 - cos(norm))*(invnorm*invnorm)*sk
				  + ((norm - sin(norm))*sk*sk)*(invnorm*invnorm*invnorm);

	  }


}



#endif /* ORB_SLAM2_THIRDPARTY_G2O_G2O_TYPES_SI3_OPS_H_ */
