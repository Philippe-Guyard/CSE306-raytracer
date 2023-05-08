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
            raytracer/bvh.hpp
            raytracer/geometry.h 
            raytracer/ray.h 
            raytracer/vector.h 
            raytracer/sphere.h
            raytracer/vector.cpp 
            raytracer/sphere.cpp
            raytracer/scene.h
            raytracer/triangle.hpp
            raytracer/mesh.hpp
            raytracer/texture.hpp)

add_executable(sphere ${SOURCES})