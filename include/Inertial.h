/*
 * Inertial.h
 *
 *  Created on: 22Nov.,2018
 *      Author: peter
 */

#ifndef ORB_SLAM2_INCLUDE_INERTIAL_H_
#define ORB_SLAM2_INCLUDE_INERTIAL_H_
#include <Eigen/Dense>

namespace ORB_SLAM2
{
struct  ImuData{
	Eigen::Vector3d accel;
	Eigen::Vector3d gyro;
	double time;
};



} // namespace ORBSLAM_2



#endif /* ORB_SLAM2_INCLUDE_INERTIAL_H_ */
