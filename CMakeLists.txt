cmake_minimum_required(VERSION 3.1)

project(gt_traj_arbitration)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

find_package(catkin REQUIRED COMPONENTS
  urdf
  cnr_controller_interface
  state_space_filters
  cnr_hardware_interface
  pluginlib
  sensor_msgs
  std_msgs
  name_sorting
  joint_limits_interface
  rospy
  std_msgs
  message_generation
  tf_conversions
)

find_package(Eigen3 3.3 REQUIRED NO_MODULE)

#add_service_files(
#  FILES
#  ioc.srv
#)

#generate_messages(
#  DEPENDENCIES
#  std_msgs
#)

catkin_package(
 INCLUDE_DIRS include
 LIBRARIES gt_traj_arbitration
 CATKIN_DEPENDS urdf cnr_controller_interface state_space_filters cnr_hardware_interface pluginlib sensor_msgs std_msgs name_sorting joint_limits_interface
 DEPENDS EIGEN3
)

include_directories( include ${catkin_INCLUDE_DIRS} )




add_library(${PROJECT_NAME} src/${PROJECT_NAME}/gt_traj_arbitration.cpp )
add_dependencies(${PROJECT_NAME} ${${PROJECT_NAME}_EXPORTED_TARGETS} ${catkin_EXPORTED_TARGETS})
target_link_libraries(${PROJECT_NAME} ${catkin_LIBRARIES} Eigen3::Eigen)
target_compile_options(${PROJECT_NAME} PUBLIC -Wall -faligned-new
   $<$<CONFIG:Release>:-Ofast -funroll-loops -ffast-math >)

if(${CMAKE_VERSION} VERSION_GREATER  "3.16.0")
target_precompile_headers(${PROJECT_NAME} PUBLIC
<cmath>
<Eigen/Core>
<cnr_logger/cnr_logger.h>
<ros/time.h>
<sensor_msgs/JointState.h>
)
endif()

add_executable(main src/gt_traj_arbitration/main.cpp)
target_link_libraries(main ${catkin_LIBRARIES})
add_dependencies(main gt_traj_arbitration)


