cmake_minimum_required(VERSION 3.10)
find_package(/home/sisyphus/.local/lib/python3.10/site-packages/pybind11/share/cmake/pybind11/pybind11 REQUIRED)
project(fast_iter)
pybind11_add_module(iterr iterr.cpp)
