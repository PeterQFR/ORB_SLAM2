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

#ifndef G2O_SI3QUAT_H_
#define G2O_SI3QUAT_H_



#include <Eigen/Core>
#include <Eigen/Geometry>

namespace g2o {
  using namespace Eigen;


  typedef Matrix<double, 6, 1> Vector6d;
  typedef Matrix<double, 16> Vector16d; //for maxvector
  //qw, qx, qy, qz, px, py, pz, vx, vy, vz, b
  typedef Matrix<double, 15> Vector15d; //For minimal vector

/** This class holds velocity information as well as pose information.
 *It is priamrily used for holding the preintegral values for imu integration.
 */
  class SI3Quat {
    public:
      EIGEN_MAKE_ALIGNED_OPERATOR_NEW;


    protected:

      SE3Quat _t;
      Vector3d _v; //velocity estimate world frame
      Vector3d _bw;
      Vector3d _ba;
    public:
      SI3Quat(){
        _v.setZero();

      }

      SI3Quat(const Quaterniond q, const Vector3d p, const Vector3d v,
    		  const Vector3d ba, const Vector3d bw):
      _t(q, p), _v(v), _ba(ba), _bw(bw){};

      inline void setBiases(Vector3d bw, Vector3d ba)
      {
    	_bw = bw;
    	_ba = ba;
      }

      inline SI3Quat operator* (const SI3Quat& tr2) const{
    	  SI3Quat result(*this);
    	  result._t *= tr2._t;
    	  result._v *= tr2.v;
    	  return result;
      }

      inline SI3Quat& operator*= (const SI3Quat& tr2){
    	  _t*=tr2._t;
    	  _v*=tr2._v;

    	  return *this;
      }
      inline Vector

      /**
       * \brief increment with the
       */
      inline SI3Quat operator* (const Vector6d& aw, const double& dt) const {
    	  Vector3d a = ag.block(0,0,1,3);
    	  Vector3d w = ag.block(0,3,1,3);
    	  SI3Quat deltaParts = preIntegrate(a, w, dt);
    	  *this*=deltaParts;
    	  return *this;
      }

      inline SI3Quat inverse() const{
    	  SI3Quat ret;
    	  ret._t=_t.inverse();
    	  ret._v=-_v;
    	  return ret;
      }


      /**
       * \brief this method determines the preintegral of SE2Quat plus velocity
       * and returns this as a SI3Quat where only the
       */
      inline void preIntegrate(const Vector3d& a, const Vector3d& w,
    		  const double& dt){

    	  Vector3d bw = this->_bw;
    	  Vector3d ba = this->_ba;
    	  Vector3d g(0.0, 0.0, -9.8);
    	  Quaterniond deltaQ= zeroOrderIntegrator((w - bw)*dt)*Quaterniond::Identity();
    	  Vector3d accelIntegral = this->_t._r.toRotationMatrix * (a-ba)*dt;
    	  Vector3d gravityIntegral = g*dt;
    	  Vector3d deltaV = accelIntegral + gravityIntegral;

    	  Vector3d deltaP = deltaV*dt + 0.5*dt*(accelIntegral + gravityIntegral);

    	  this->_t *=SE3Quat(deltaQ, deltaP);
    	  this->_v += deltaV;

      }
      inline const Vector3d& velocity() const {return _v;}

      inline void setVelocity(const Vector3d& v_) {_v=v_;}



  };

  inline std::ostream& operator <<(std::ostream& out_str, const SE3Quat& se3)
  {
    out_str << se3.to_homogeneous_matrix()  << std::endl;
    return out_str;
  }

} // end namespace

#endif
