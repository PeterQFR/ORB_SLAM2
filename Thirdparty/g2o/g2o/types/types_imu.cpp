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


} //ns g2o

