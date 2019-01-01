/*
 * types_imu.cpp
 *
 *  Created on: 11Nov.,2018
 *      Author: peter
 */


#include "types_imu.h"
#include "../core/factory.h"

namespace g2o {

using namespace std;
using namespace Eigen;

G2O_REGISTER_TYPE_GROUP(types_imu);
G2O_REGISTER_TYPE(VERTEX_VELOCITY, VertexVelocity);
G2O_REGISTER_TYPE(VERTEX_BIAS, VertexBias);
G2O_REGISTER_TYPE(EDGE_IMU_UPDATE, EdgeImuUpdate);
G2O_REGISTER_TYPE(EDGE_BIAS_UPDATE, EdgeBiasUpdate);

bool VertexVelocity::read(std::istream& is){
	Vector3d  velocity;
	for (int i=0; i< 3; i++)
	{
		is >> velocity[i];
	}
	setEstimate(velocity);
}

bool VertexVelocity::write(std::ostream& os) const {
	Vector3d velocity = estimate();
	for (int i = 0; i < 3; i++)
	{
		os << velocity[i] << " ";
	}
	return os.good();
}

bool VertexBias::read(std::istream& is){
	Vector6d  velocity;
	for (int i=0; i< 6; i++)
	{
		is >> velocity[i];
	}
	setEstimate(velocity);
}

bool VertexBias::write(std::ostream& os) const {
	Vector6d velocity = estimate();
	for (int i = 0; i < 6; i++)
	{
		os << velocity[i] << " ";
	}
	return os.good();
}


bool EdgeBiasUpdate::read(std::istream& is){
	for (int i=0; i< 6; i++)
	{
		is >> _measurement[i];
	}
	for (int i=0; i<6; i++)
	{
		for (int j=i; j<6; j++)
		{
			is >> information()(i,j);
			if (i!=j)
				information()(j,i)=information()(i,j);
		}
	}
	return true;
}

bool EdgeBiasUpdate::write(std::ostream& os) const{
	 for (int i=0; i<6; i++){
	    os << measurement()[i] << " ";
	  }

	  for (int i=0; i<6; i++)
	    for (int j=i; j<6; j++){
	      os << " " <<  information()(i,j);
	    }
	  return os.good();
	return true;
}

bool EdgeImuUpdate::read(std::istream& is){

	return false;
}

bool EdgeImuUpdate::write(std::ostream& os) const{
	return false;
}
} //ns g2o

