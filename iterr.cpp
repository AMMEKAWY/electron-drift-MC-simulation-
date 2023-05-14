#include <iostream>
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
//#include "xtensor-python/pyvectorize.hpp"
#include <vector>
#include <pybind11/stl.h>


namespace py = pybind11;

//py::array_t<double, py::array::f_style> arr({3});


py::array_t<double, py::array::c_style> arr({3});
double* data;

std::vector<double> iterr(double &x, double &y, double &z, double &xd, double &yd, double &zd, double &xd2, double &yd2, double &zd2, double &vx, double &vy, double &vz, double &delt, double &h){

	//static double arr[3]={0,0,0};
	//py::list arr[3]={0,0,0};
	
	//auto arr= py::array_t<double>(3);
	
	std::vector<double> arr;
	
	while (x < xd){
	
	
		vx=vx+delt*h;
		x=vx*h+delt*h*h/2+x;
		y=y+vy*h;
		z=z+vz*h;
		
		
		
		if (x > xd2 || y < yd || y > yd2 || z < zd || z > zd2 ) {
		
			break;
		
		}
		
			
	}
	
	//double arr[3]={x, y, z};
	
	//py::list::append arr=x;
	
	//py::array_t<double, py::array::c_style> arr({});
		
	//arr[3]={1,2,3};
	
	//arr.data(2);
	
	//auto data = arr.data(x, y, z);	
	
	//py::array_t<double, py::array::c_style> arr({});
	
	//py::array::c_style> arr({3}, data);
	
	//arr=({}, data)
	
	//py::array::c_style> arr({3}, data);

	arr={x, y, z, vx, vy, vz};
	
	return arr;
		
}

PYBIND11_MODULE(iterr, handle){

	

	handle.def("iterr",&iterr);

	
}
