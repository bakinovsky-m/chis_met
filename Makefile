SOURCES = main.cpp header.cpp
CXX_FLAGS += -g -Wall -Wextra -Wpedantic

all: compile
	./build/main

compile:
	mkdir -p build
	g++ $(SOURCES) -o build/main $(CXX_FLAGS)

clean:
	rm -rf build