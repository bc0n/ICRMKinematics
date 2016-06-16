# ICRMKinematics
Solving the forward and inverse kinematics of the 5D maipulator.
The main project _ICRMKinematics_ includes methods for searching the joint and kinematic parameter space, encapsulating these into a dll. Sub-projects test various aspects of this dll.

## Notes on configuring visual studio
### ICRMKinematics
This builds an x86 dll callable by the outer LabVIEW loop. As such, the project properties are:
* Configuration: All Configurations, Platform: Win32 (should be Active)
* General -- Configuration Type = dll
* C/C++ -- All Options -- Additional Include Directories = D:\flat\interleavedCatheter\eigen-eigen-bdd17ee3b1b3\; D:\flat\interleavedCatheter\nlopt-2.4.2-dll32;
* C/C++ -- All Options -- Enable Enhanced Instruction Set = SIMD 2
* C/C++ -- All Options -- Open MP Support = Yes
* Linker -- All Options -- Additional Dependencies = libnlopt-0.lib
* Linker -- All Options -- Additional Library Directories = D:\flat\interleavedCatheter\nlopt-2.4.2-dll32;
