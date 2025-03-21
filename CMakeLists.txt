cmake_minimum_required(VERSION 3.10)
project(ga-ddp)

set(CMAKE_BUILD_TYPE "Release")
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -Wall -g")

# find dependencies
find_package(ament_cmake REQUIRED)
find_package(rclcpp REQUIRED)
find_package(rclpy REQUIRED)
find_package(Eigen3 REQUIRED)
find_package(OpenMP REQUIRED)
find_package(sensor_msgs REQUIRED)
find_package(trajectory_msgs REQUIRED)

include_directories(
    include
    ${EIGEN3_INCLUDE_DIR}
)

add_library(GA
    src/GA/algebra.cpp
    src/GA/robotics.cpp
)

add_library(OC
    src/OC/ilqg.cpp
    src/OC/boxqp.cpp
)
target_link_libraries(OC OpenMP::OpenMP_CXX)

add_executable(GA_test test/GA_test.cpp)
target_link_libraries(GA_test GA)

add_executable(iLQG_linear_test test/iLQG_linear_test.cpp)
target_link_libraries(iLQG_linear_test GA OC OpenMP::OpenMP_CXX)

add_executable(iLQG_robot_test test/iLQG_robot_test.cpp)
target_link_libraries(iLQG_robot_test GA OC OpenMP::OpenMP_CXX)

# add_executable(robot_mpc_test tests/robot_mpc_test.cpp)
# ament_target_dependencies(robot_mpc_test rclcpp sensor_msgs trajectory_msgs)
# target_link_libraries(robot_mpc_test ${COMMON_LIBS})

# Install targets
install(TARGETS
    GA_test
    iLQG_linear_test
    iLQG_robot_test
#     robot_mpc_test
  DESTINATION lib/${PROJECT_NAME}
)

if(BUILD_TESTING)
    find_package(ament_lint_auto REQUIRED)
    set(ament_cmake_cpplint_FOUND TRUE)
    ament_lint_auto_find_test_dependencies()
endif()

ament_package()