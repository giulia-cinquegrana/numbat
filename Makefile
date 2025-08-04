OMP_CXX     = /usr/local/opt/llvm/bin/clang++
OMP_FLAGS   = -std=c++17 -O2 -fopenmp -I/usr/local/opt/libomp/include
OMP_LDFLAGS = -L/usr/local/opt/libomp/lib

SRC         = body.cpp
OUT         = body

all: $(OUT)

$(OUT): $(SRC)
	$(OMP_CXX) $(OMP_FLAGS) $(SRC) -o $(OUT) $(OMP_LDFLAGS)

clean:
	rm -f $(OUT)
