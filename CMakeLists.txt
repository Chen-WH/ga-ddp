cmake_minimum_required(VERSION 3.10)
project(ga-ddp)

# 设置 C++ 标准
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# 添加编译选项
add_compile_options(-O3 -fopenmp)

# find dependencies
find_package(Eigen3 REQUIRED)
find_package(rclcpp REQUIRED)
find_package(sensor_msgs REQUIRED)
find_package(trajectory_msgs REQUIRED)
find_package(ament_cmake REQUIRED)
find_package(OpenMP REQUIRED)

# 设置头文件路径
include_directories(${EIGEN3_INCLUDE_DIR})
include_directories(${CMAKE_SOURCE_DIR}/include)

# 添加源代码
add_library(GA STATIC
    src/GA_algebra.cpp
    src/GA_robot.cpp
)

add_library(OC STATIC
    src/ilqr.cpp
    src/boxqp.cpp
)

# 设置库的输出路径
set(LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR}/lib)

# 设置公共链接库
set(COMMON_LIBS GA OC OpenMP::OpenMP_CXX)

# 安装目录
install(DIRECTORY
  #launch
  DESTINATION share/${PROJECT_NAME}/
)

# 示例应用
add_executable(GA_algebra_test tests/GA_algebra_test.cpp)
target_link_libraries(GA_algebra_test ${COMMON_LIBS})

add_executable(GA_robot_test tests/GA_robot_test.cpp)
target_link_libraries(GA_robot_test ${COMMON_LIBS})

add_executable(iLQR_boxqp_test tests/iLQR_boxqp_test.cpp)
target_link_libraries(iLQR_boxqp_test ${COMMON_LIBS})

add_executable(iLQR_linear_test tests/iLQR_linear_test.cpp)
target_link_libraries(iLQR_linear_test ${COMMON_LIBS})

add_executable(iLQR_robot_test tests/iLQR_robot_test.cpp)
target_link_libraries(iLQR_robot_test ${COMMON_LIBS})

add_executable(robot_mpc_test tests/robot_mpc_test.cpp)
ament_target_dependencies(robot_mpc_test rclcpp sensor_msgs trajectory_msgs)
target_link_libraries(robot_mpc_test ${COMMON_LIBS})

# 安装可执行文件
install(TARGETS
    GA_algebra_test
    GA_robot_test
    iLQR_boxqp_test
    iLQR_linear_test
    iLQR_robot_test
    robot_mpc_test
  DESTINATION lib/${PROJECT_NAME}
)

if(BUILD_TESTING)
  find_package(ament_lint_auto REQUIRED)
  set(ament_cmake_cpplint_FOUND TRUE)
  ament_lint_auto_find_test_dependencies()
endif()

ament_package()