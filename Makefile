CXXFLAGS = -std=c++0x -O3

heatmap: heatmap.o
	g++ $(LDFLAGS) $^ -o $@

clean:
	rm heatmap heatmap.o
