
#include <fstream>
#include <iostream>
#include "kinematicsDLL.h"

void main(){
	//load test file
	//char fname[] = "testTriangleData.dat"; // n x 38 double 
	//char fname[] = "binary_big.dat";
	char fname[] = "binary_little.dat";
	int ncols = 4;//38;
	int nrows = 0;
	long fileSize = 0;
	char *buffer;
	double rowData[38];//ncols

	//open
	printf("\n\nOpening %s\n", fname);
	std::fstream fs;
	fs.open(fname, std::fstream::in | std::fstream::binary); //name, read/write mode|binary
	if (!fs) { std::cerr << "warning, could not open " << fname << std::endl; return; };
	//nrows
	fs.seekg(0,fs.end);
	fileSize = fs.tellg();
	nrows = fileSize / sizeof(double)/ncols;
	fs.seekg(0, fs.beg);
	printf("fileSize %d with %d rows\n", fileSize, nrows);
	//assign
	buffer = new char[fileSize];
	fs.read(buffer, fileSize);
	printf("read %d chars\n", fs.gcount());
	fs.close();

	//interp
	double *pd = reinterpret_cast<double*>(buffer);
	double rowData[40];
	for (int i = 0; i < ncols*nrows; i++) {
		rowData[i] = *(pd + i);
	}
	printf("%f %f %f %f\n", rowData[0], rowData[1], rowData[2], rowData[3]);
	printf("%f %f %f %f\n", rowData[4], rowData[5], rowData[6], rowData[7]);
	
	delete[] buffer; //everybody do your share
	
	
	std::cout << "\n\nPress Enter";
	std::getchar();
	printf("done\n");

}

