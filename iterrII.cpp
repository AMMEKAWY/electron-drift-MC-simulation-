#include <iostream>
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
//#include "xtensor-python/pyvectorize.hpp"
#include <vector>
#include <pybind11/stl.h>
#include <math.h>

namespace py = pybind11;

//py::array_t<long double, py::array::f_style> arr({3});


py::array_t<long double, py::array::c_style> arr({3});
long double* data;

std::vector<long double> iterr(long double x, long double y, long double z, long double xd, long double yd, long double zd, long double xd2, long double yd2, long double zd2, long double vx, long double vy, long double vz, long double delt, long double h){

	//static long double arr[3]={0,0,0};
	//py::list arr[3]={0,0,0};
	
	//auto arr= py::array_t<long double>(3);
	
	std::vector<long double> arr;
	
	double xr=x;
	double yr=y;
	double zr=z;
	
	double r=sqrt((x-xr)*(x-xr)+(y-yr)*(y-yr)+(z-zr)*(z-zr));
	double rd=sqrt(xd*xd+yd*yd+zd*zd);
		
	double rat=r/rd;	
	
	while (rat<1){
	
	
		vx=vx+delt*h;
		x=vx*h+delt*h*h/2+x;
		y=y+vy*h;
		z=z+vz*h;
		
		//std::cout<<x<<std::endl;	
		
		r=sqrt((x-xr)*(x-xr)+(y-yr)*(y-yr)+(z-zr)*(z-zr));
		rd=sqrt(xd*xd+yd*yd+zd*zd);
		
		rat=r/rd;	
			
	}
	
	//long double arr[3]={x, y, z};
	
	//py::list::append arr=x;
	
	//py::array_t<long double, py::array::c_style> arr({});
		
	//arr[3]={1,2,3};
	
	//arr.data(2);
	
	//auto data = arr.data(x, y, z);	
	
	//py::array_t<long double, py::array::c_style> arr({});
	
	//py::array::c_style> arr({3}, data);
	
	//arr=({}, data)
	
	//py::array::c_style> arr({3}, data);

	arr={x, y, z, vx, vy, vz};
	
	
	
	return arr;
		
}

PYBIND11_MODULE(iterr, handle){

	

	handle.def("iterr",iterr);

	
}
