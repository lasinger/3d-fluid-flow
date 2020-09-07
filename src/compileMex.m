%%% compile mex files

% -D__QUIET__ 
% no openmp: -D_NO_OPENMP
% OPTIMFLAGS = /O2 /Oy- /DNDEBUG /openmp
% for more speed compile externally with other compiler

% download Eigen library first:

% with openmp:
windows=1;
if windows
    disp('windows');
     mex -largeArrayDims CC=g++ CXX=g++ LD=g++ -v OPTIMFLAGS="$OPTIMFLAGS /O2 /DNDEBUG /openmp" CXXFLAGS="$CXXFLAGS -fPIC -fopenmp -march=native" LDFLAGS="$LDFLAGS -fopenmp" -O -I./ -output ./render2dallcammex ./Source/Render2AllCam_Mex.cpp
      mex -largeArrayDims CC=g++ CXX=g++ LD=g++ -v OPTIMFLAGS="$OPTIMFLAGS /O2 /DNDEBUG /openmp" CXXFLAGS="$CXXFLAGS -fPIC -fopenmp -march=native" LDFLAGS="$LDFLAGS -fopenmp" -O -I./ -output ./gradientmex ./Source/Gradient_Mex.cpp
      mex -largeArrayDims CC=g++ CXX=g++ LD=g++ -v OPTIMFLAGS="$OPTIMFLAGS /O2 /DNDEBUG /openmp" CXXFLAGS="$CXXFLAGS -fPIC -fopenmp -march=native" LDFLAGS="$LDFLAGS -fopenmp" -O -I./ -output ./gradientpolymex ./Source/GradPoly_Mex.cpp
      mex -largeArrayDims CC=g++ CXX=g++ LD=g++ -v OPTIMFLAGS="$OPTIMFLAGS /O2 /DNDEBUG /openmp" CXXFLAGS="$CXXFLAGS -fPIC -fopenmp -march=native" LDFLAGS="$LDFLAGS -fopenmp" -O -I./ -output ./gradientpolyPartmex ./Source/GradPolyPart_Mex.cpp
      mex -largeArrayDims CC=g++ CXX=g++ LD=g++ -v OPTIMFLAGS="$OPTIMFLAGS /O2 /DNDEBUG /openmp" CXXFLAGS="$CXXFLAGS -fPIC -fopenmp -march=native" LDFLAGS="$LDFLAGS -fopenmp" -O -I./ -output ./gradientintmex ./Source/GradInt_Mex.cpp
      mex -largeArrayDims CC=g++ CXX=g++ LD=g++ -v OPTIMFLAGS="$OPTIMFLAGS /O2 /DNDEBUG /openmp" CXXFLAGS="$CXXFLAGS -fPIC -fopenmp -march=native" LDFLAGS="$LDFLAGS -fopenmp" -O -I./ -output ./interpFlowToPartmex ./Source/InterpFlowToPart_Mex.cpp
      mex -largeArrayDims CC=g++ CXX=g++ LD=g++ -v OPTIMFLAGS="$OPTIMFLAGS /O2 /DNDEBUG /openmp" CXXFLAGS="$CXXFLAGS -fPIC -fopenmp -march=native" LDFLAGS="$LDFLAGS -fopenmp" -O -I./ -output ./derivativeFlowPartmex ./Source/DerivativeFlowPart_Mex.cpp
      mex -largeArrayDims CC=g++ CXX=g++ LD=g++ -v OPTIMFLAGS="$OPTIMFLAGS /O2 /DNDEBUG /openmp" CXXFLAGS="$CXXFLAGS -fPIC -fopenmp -march=native" LDFLAGS="$LDFLAGS -fopenmp" -O -I./ -output ./derivativePart2Gridmex ./Source/DerivativePart2Grid_Mex.cpp
      mex -largeArrayDims CC=g++ CXX=g++ LD=g++ -v OPTIMFLAGS="$OPTIMFLAGS /O2 /DNDEBUG /openmp" CXXFLAGS="$CXXFLAGS -fPIC -fopenmp -march=native" LDFLAGS="$LDFLAGS -fopenmp" -O -I./ -I./Source/eigen/ -output ./triangulatePartmex ./Source/Triang_Mex.cpp
      mex -largeArrayDims CC=g++ CXX=g++ LD=g++ -v OPTIMFLAGS="$OPTIMFLAGS /O2 /DNDEBUG /openmp" CXXFLAGS="$CXXFLAGS -fPIC -fopenmp -march=native" LDFLAGS="$LDFLAGS -fopenmp" -O -I./ -I./Source/eigen/ -output ./proj2dmex ./Source/Proj2d_Mex.cpp
      mex -largeArrayDims CC=g++ CXX=g++ LD=g++ -v OPTIMFLAGS="$OPTIMFLAGS /O2 /DNDEBUG /openmp" CXXFLAGS="$CXXFLAGS -fPIC -fopenmp -march=native" LDFLAGS="$LDFLAGS -fopenmp" -O -I./ -I./Source/eigen/ -output ./triangulatePartPolymex ./Source/TriangPoly_Mex.cpp
else
    disp('linux');
    % -march=native maybe not so wise on Euler (different cpus)
    mex -largeArrayDims CC=g++ CXX=g++ LD=g++ -v OPTIMFLAGS='$OPTIMFLAGS /O2 /DNDEBUG /openmp' CXXFLAGS='-fPIC -fopenmp -O2' LDFLAGS='\$LDFLAGS -fopenmp' -lgomp -I./ -output ./render2dallcammex ./Source/Render2AllCam_Mex.cpp
    mex -largeArrayDims CC=g++ CXX=g++ LD=g++ -v OPTIMFLAGS='$OPTIMFLAGS /O2 /DNDEBUG /openmp' CXXFLAGS='-fPIC -fopenmp -O2' LDFLAGS='\$LDFLAGS -fopenmp' -lgomp -I./ -output ./gradientmex ./Source/Gradient_Mex.cpp
    mex -largeArrayDims CC=g++ CXX=g++ LD=g++ -v OPTIMFLAGS='$OPTIMFLAGS /O2 /DNDEBUG /openmp' CXXFLAGS='-fPIC -fopenmp -O2' LDFLAGS='\$LDFLAGS -fopenmp' -lgomp -I./ -output ./gradientpolymex ./Source/GradPoly_Mex.cpp
    mex -largeArrayDims CC=g++ CXX=g++ LD=g++ -v OPTIMFLAGS='$OPTIMFLAGS /O2 /DNDEBUG /openmp' CXXFLAGS='-fPIC -fopenmp -O2' LDFLAGS='\$LDFLAGS -fopenmp' -lgomp -I./ -output ./gradientpolyPartmex ./Source/GradPolyPart_Mex.cpp
    mex -largeArrayDims CC=g++ CXX=g++ LD=g++ -v OPTIMFLAGS='$OPTIMFLAGS /O2 /DNDEBUG /openmp' CXXFLAGS='-fPIC -fopenmp -O2' LDFLAGS='\$LDFLAGS -fopenmp' -lgomp -I./ -output ./gradientintmex ./Source/GradInt_Mex.cpp
    mex -largeArrayDims CC=g++ CXX=g++ LD=g++ -v OPTIMFLAGS='$OPTIMFLAGS /O2 /DNDEBUG /openmp' CXXFLAGS='-fPIC -fopenmp -O2' LDFLAGS='\$LDFLAGS -fopenmp' -lgomp -I./ -output ./interpFlowToPartmex ./Source/InterpFlowToPart_Mex.cpp
    mex -largeArrayDims CC=g++ CXX=g++ LD=g++ -v OPTIMFLAGS='$OPTIMFLAGS /O2 /DNDEBUG /openmp' CXXFLAGS='-fPIC -fopenmp -O2' LDFLAGS='\$LDFLAGS -fopenmp' -lgomp -I./ -output ./derivativeFlowPartmex ./Source/DerivativeFlowPart_Mex.cpp
    mex -largeArrayDims CC=g++ CXX=g++ LD=g++ -v OPTIMFLAGS='$OPTIMFLAGS /O2 /DNDEBUG /openmp' CXXFLAGS='-fPIC -fopenmp -O2' LDFLAGS='\$LDFLAGS -fopenmp' -lgomp -I./ -output ./derivativePart2Gridmex ./Source/DerivativePart2Grid_Mex.cpp
    mex -largeArrayDims CC=g++ CXX=g++ LD=g++ -v OPTIMFLAGS='$OPTIMFLAGS /O2 /DNDEBUG /openmp' CXXFLAGS='-fPIC -fopenmp -O2' LDFLAGS='\$LDFLAGS -fopenmp' -lgomp -I./ -I./Source/eigen/ -output ./triangulatePartmex ./Source/Triang_Mex.cpp
    mex -largeArrayDims CC=g++ CXX=g++ LD=g++ -v OPTIMFLAGS='$OPTIMFLAGS /O2 /DNDEBUG /openmp' CXXFLAGS='-fPIC -fopenmp -O2' LDFLAGS='\$LDFLAGS -fopenmp' -lgomp -I./ -I./Source/eigen/ -output ./proj2dmex ./Source/Proj2d_Mex.cpp
    mex -largeArrayDims CC=g++ CXX=g++ LD=g++ -v OPTIMFLAGS='$OPTIMFLAGS /O2 /DNDEBUG /openmp' CXXFLAGS='-fPIC -fopenmp -O2' LDFLAGS='\$LDFLAGS -fopenmp' -lgomp -I./ -I./Source/eigen/ -output ./triangulatePartPolymex ./Source/TriangPoly_Mex.cpp
end
