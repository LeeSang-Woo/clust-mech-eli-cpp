# clust-mech-eli-cpp
 This C++ code was used in the paper "Two types of critical cell density for mechanical elimination of abnormal cell clusters from epithelial tissue" for the simulations of vertex dynamics model.

### Environment
- OS : Windows (tested with Windows 10)
- Compiler : Visual C++ (tested with Visual C++ 2022)
- IDE : Visual Studio (tested with Visual Studio Community 2022 64-bit)

### Build instruction (Windows)
1. Download Visual Studio Community (https://visualstudio.microsoft.com/).
2. Install Visual Studio Community. Check "Desktop development with C++", this will install Visual C++ compiler.
3. Create the folder and download the source codes.
4. Launch Visual Studio and open the solution file (Menu bar> File> Open> Project/Solution> Select "clust-mech-eli.sln").
  
![vs01](https://user-images.githubusercontent.com/46122090/156304341-9d2112b5-8999-45af-8f90-1c3d2d2f633d.png)

5. Change Solution Configrations to "Release".

![vs03](https://user-images.githubusercontent.com/46122090/156312753-5fcfd191-5506-4ba1-a9be-479d00abb594.png)

6. Build the solution (Menu bar> Build> Build Solution).

![vs02](https://user-images.githubusercontent.com/46122090/156305191-765645de-e7c5-45a9-9c1a-16a47198c74f.png)

7. Executive file "clust-mech-eli.exe" will be built into "downloaded folder/x64/Release/".

### Execute simulation (Windows)
- For the simulation, the executive file "clust-mech-eli.exe" needs two files describing the initial condition ("initial2D_RD05_in1000_c_L####G####.csv" and "initial2D_RD05_in1000_p_L####G####.csv") and three files describing the simulation conditions ("DCModel_container_info.txt", "DCModel_param_abnormal.txt", and "DCModel_param_basic.txt").
- All these files must be in the same folder.
- When you run the simulation, "data2D_N####mu####.csv" files will be created as results.
