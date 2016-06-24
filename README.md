# ICRMKinematics
This project solves the forward and inverse kinematics of the 5D maipulator.
The main project _ICRMKinematics\ICRMKinematics_ includes methods for searching the joint and kinematic parameter space, encapsulating these into a dll.
Sub-projects test various aspects of this dll in c++, labview, and matlab.

# Configuration notes:
This project depends on [Eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page) and the [nonlinear optimization library](https://github.com/stevengj/nlopt) and assumes these are placed in the same parent folder as ICRMKinematics.
You will need to create the nlopt .lib from the included .def; see the nlopt windows readme.

## ICRMKinematics\ICRMKinematics
This contains the forward and inverse kinematics routines. It builds an x86 dll callable by the outer LabVIEW loop. As such, the project properties are:
* Configuration: All Configurations, Platform: Win32 (should be Active)
* General -- Configuration Type = dll
* C/C++ -- All Options -- Additional Include Directories = ..\..\eigen-eigen-bdd17ee3b1b3\; ..\..\nlopt-2.4.2-dll32;
* C/C++ -- All Options -- Enable Enhanced Instruction Set = SIMD 2
* C/C++ -- All Options -- Open MP Support = Yes
* Linker -- All Options -- Additional Dependencies = libnlopt-0.lib
* Linker -- All Options -- Additional Library Directories = ..\..\nlopt-2.4.2-dll32;

## ICRMKinematics\ICRMKinematics_testDLL
Simple tests of the functions in the exported dll.

## ICRMKinematics\matlab
Contains a matlab implementation of the various methods as well as functions for examining the resulting data.
### ICRMKinematics\matlab\testDLLx64.m
For testing the ICRMKinematics.dll from matlab. My setup requires building an x64 ICRMKinematics.dll, which then needs the x64 nlopt.  For running, the nloptx64 needs to be in the same folder as ICRMKinematics.dll x64.