CXXFLAGS = -std=c++0x -O3 $(shell pkg-config --cflags Magick++)
LDFLAGS = $(shell pkg-config --libs Magick++)

heatmap: heatmap.o
	g++ $^ -o $@ $(LDFLAGS)

clean:
	rm heatmap heatmap.o
