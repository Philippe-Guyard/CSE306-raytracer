cmake_minimum_required(VERSION 3.5)

project(sphere)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_FLAGS "-Wall -Wextra -fopenmp")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

set(SOURCES main.cpp
            mesh_reader.hpp
            raytracer/geometry.h 
            raytracer/ray.h 
            raytracer/vector.h 
            raytracer/sphere.h
            raytracer/vector.cpp 
            raytracer/sphere.cpp
            raytracer/scene.h
            raytracer/triangles.h)

add_executable(sphere ${SOURCES})