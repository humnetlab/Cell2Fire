# Cell2Fire + ML
## CP: Kept for Minho's reference

### cmake build
Under the folder : 2_BehavePlus/cell2fire_BP
- cd /Users/minho/Documents/GitHub/c3ai-fire/2_BehavePlus/cell2fire_BP/Cell2FireC (Use full path)
- mkdir build
- cd build
- cmake ..
- make

#### Debugging
(1) If you get **CMake Error at CMakeLists.txt:17 (find_package):**
# cmake -DCMAKE_PREFIX_PATH /home/minho/fires/cell2fireML/libtorch ..
# cmake 

(2) If you get **PermissionError: [Errno 13] Permission denied:**
Re-run cmake build to create cmake files again. Then try "make"

**NOTE: You need to download libtorch and place it in the folder "cell2fire_BP". The current libtorch directory is *empty*. Download libtorch at pytorch.org. Use the following specifications:**
- PyTorch Build: "Preview (Nightly)" 
- OS: Linux / Windows / Mac
- Package: "Libtorch"
- Language: C++/Java
- Compute Platform: Linux(CUDA11.6 or CPU) / Windows(Debug version!) / Mac(Default)

### Execute Cell2Fire
- mv Cell2Fire ..
# The executable file has to be moved from the "build" folder to the Cell2FireC folder
- cd ../..
- python main.py --input-instance-folder data_BP/f101/ --output-folder results/BP_f101_test/ --ignitions --sim-years 1 --nsims 1 --grids --finalGrid --weather rows --nweathers 1 --Fire-Period-Length 1.0 --ROS-CV 0.0 --output-messages --seed 123 --IgnitionRad 1 --stats --verbose --allPlots