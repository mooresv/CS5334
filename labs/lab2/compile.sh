export PFETCH=2
icc  -march=core-avx2 -qno-offload -DAOS -DDOUBLEPREC -qopenmp -O3 -g -qopt-prefetch=$PFETCH icpmain.cpp eigend.cpp dtime.c -o icpdblaosavx.exe
icc  -xSSE2 -qno-offload -DAOS -DDOUBLEPREC -qopenmp -O3 -g -qopt-prefetch=$PFETCH icpmain.cpp eigend.cpp dtime.c -o icpdblaossse.exe
icc  -march=core-avx2 -qno-offload -DDOUBLEPREC -qopenmp -O3 -g -qopt-prefetch=$PFETCH icpmain.cpp eigend.cpp dtime.c -o icpdblsoaavx.exe
icc  -xSSE2 -qno-offload -DDOUBLEPREC -qopenmp -O3 -g -qopt-prefetch=$PFETCH icpmain.cpp eigend.cpp dtime.c -o icpdblsoasse.exe
icc  -march=core-avx2 -qno-offload -DAOS -qopenmp -O3 -g -qopt-prefetch=$PFETCH icpmain.cpp eigend.cpp dtime.c -o icpfltaosavx.exe
icc  -xSSE2 -qno-offload -DAOS -qopenmp -O3 -g -qopt-prefetch=$PFETCH icpmain.cpp eigend.cpp dtime.c -o icpfltaossse.exe
icc  -march=core-avx2 -qno-offload -qopenmp -O3 -g -qopt-prefetch=$PFETCH icpmain.cpp eigend.cpp dtime.c -o icpfltsoaavx.exe
icc  -xSSE2 -qno-offload -qopenmp -O3 -g -qopt-prefetch=$PFETCH icpmain.cpp eigend.cpp dtime.c -o icpfltsoasse.exe
icc  -no-vec -no-fma -qno-offload -DAOS -DDOUBLEPREC -qopenmp -O3 -g -qopt-prefetch=$PFETCH icpmain.cpp eigend.cpp dtime.c -o icpdblaos.exe
icc  -no-vec -no-fma -qno-offload -DDOUBLEPREC -qopenmp -O3 -g -qopt-prefetch=$PFETCH icpmain.cpp eigend.cpp dtime.c -o icpdblsoa.exe
icc  -no-vec -no-fma -qno-offload -DAOS -qopenmp -O3 -g -qopt-prefetch=$PFETCH icpmain.cpp eigend.cpp dtime.c -o icpfltaos.exe
icc  -no-vec -no-fma -qno-offload -qopenmp -O3 -g -qopt-prefetch=$PFETCH icpmain.cpp eigend.cpp dtime.c -o icpfltsoa.exe
