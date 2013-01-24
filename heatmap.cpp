#include <cassert>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <iostream>
#include <fstream>
#include <functional>
#include <map>
#include <set>
#include <vector>
#include <Magick++.h>

struct Point {
    Point(double lat, double lon, bool radian=true) {
        if (radian) {
            lat_r = lat;
            lon_r = lon;
            lat_d = 180. * lat /  M_PI;
            lon_d = 180. * lon /  M_PI;
        }
        else {
            lat_d = lat;
            lon_d = lon;
            lat_r = M_PI * lat /  180.;
            lon_r = M_PI * lon /  180.;
        }
    }

    double lat_r, lon_r;
    double lat_d, lon_d;
};

class TrackRaster {
    public:
        TrackRaster() { }
        void mark(int x, int y) { data.insert(std::pair<int, int>(x, y)); }
        bool is_marked(int x, int y) { return data.find(std::pair<int, int>(x,y)) != data.end(); }
        std::set<std::pair<int, int>>::const_iterator begin() const { return data.begin(); }
        std::set<std::pair<int, int>>::const_iterator end() const { return data.end(); }
    private:
        std::set<std::pair<int, int>> data;
};

TrackRaster make_grid(const std::vector<Point>& track, std::function<std::pair<int, int>(const Point&)> conv, int radius) {
    TrackRaster tr;
    std::pair<int,int> prev(-1, -1);
    for (auto p : track) {
        auto coords = conv(p);
        if (prev.first >= 0) {
            // Bresenham's line algorithm
            int x0 = prev.first;
            int y0 = prev.second;
            int x1 = coords.first;
            int y1 = coords.second;
            int dx = abs(x1 - x0);
            int dy = abs(y1 - y0);
            int sx = (x0 < x1) ? 1 : -1;
            int sy = (y0 < y1) ? 1 : -1;
            int err = dx - dy;

            while (true) {
                tr.mark(x0, y0);
                if (x0 == x1 && y0 == y1)
                    break;
                int e2 = 2 * err;
                if (e2 > -dy) {
                    err -= dy;
                    x0 += sx;
                }
                if (e2 < dx) {
                    err += dx;
                    y0 += sy;
                }
            }
        }
        prev = coords;
    }
    return std::move(tr);
}

void write_color(int d, int low, int high, unsigned char *rgb) {
    if (d > high)
        d = high;

    // Wrote the reverse of the color scheme I wanted, so just did 1 -
    // instead of recomputing...
    float v = 1 - ((float)(d - low)) / (high - low);

    if (d == 0) {
        rgb[0] = rgb[1] = rgb[2] = 0;
    }
    else if (v <= 0.25) {
        rgb[0] = 255;
        rgb[1] = 255 * (v / 0.25);
        rgb[2] = 0;
    }
    else if (v <= 0.5) {
        rgb[0] = 255 * ((0.25 - v) / 0.25);
        rgb[1] = 255;
        rgb[2] = 0;
    }
    else if (v <= 0.75) {
        rgb[0] = 0;
        rgb[1] = 255;
        rgb[2] = 255 * ((v - 0.50) / 0.25);
    }
    else {
        rgb[0] = 0;
        rgb[1] = 255 * ((0.75 - v) / 0.25);
        rgb[2] = 255;
    }
}

void draw(char *fn, const std::vector<std::vector<Point>>& tracks, const Point& p1, const Point& p2, int width, int height) {
    // dimensions of a pixel assuming grid coordinates
    double h_thresh = (fabs(p2.lat_r - p1.lat_r)) / height;
    double w_thresh = (fabs(p2.lon_r - p1.lon_r)) / width;

    // This code is soooo incorrect, but wooo approximation
    double d_lat = fabs(p2.lat_r - p1.lat_r) / height;
    double d_lon = fabs(p2.lon_r - p1.lon_r) / width;
    double max_lat = std::max(p2.lat_r, p1.lat_r);
    double min_lon = std::min(p2.lon_r, p1.lon_r);

    std::vector<TrackRaster> rasters;

    auto conv = [min_lon, max_lat, h_thresh, w_thresh] (const Point& p) -> std::pair<int, int> {
        return std::pair<int, int>((p.lon_r - min_lon) / w_thresh, (max_lat - p.lat_r) / h_thresh);
    };

    std::clog << "Processing tracks..." << std::endl;
    for (auto&& track : tracks) {
        TrackRaster raster(make_grid(track, conv, 0));
        rasters.push_back(std::move(raster));
    }

    std::clog << "Done processing" << std::endl;

    std::map<std::pair<int, int>, unsigned char> summary;
    std::clog << "Summarizing..." << std::endl;
    for (auto raster : rasters) {
        for (auto p : raster) {
            if (p.first < 0 || p.first >= width ||
                    p.second < 0 || p.second >= height)
                continue;
            ++summary[p];
        }
    }
    std::map<unsigned char, int> freq;
    for (auto ent : summary) {
        ++freq[ent.second];
    }
    // Consider 10% to be the hottest
    int cut_off = .9 * summary.size();
    int hot_val = rasters.size();
    int cold_val = freq.begin()->first;
    int elts = summary.size();
    for (auto itr = freq.rbegin(); itr != freq.rend(); ++itr) {
        elts -= itr->second;
        if (elts < cut_off) {
            hot_val = itr->first;
            break;
        }
    }

    std::clog << "Done summarizing" << std::endl;

    std::clog << "Writing..." << std::endl;
    unsigned char *pixels = new unsigned char[height * width * 3];
    unsigned char *p = pixels;
    memset(pixels, 0, height * width * 3);
    for (auto ent : summary) {
        auto x = ent.first.first;
        auto y = ent.first.second;
        assert(x < width);
        assert(y < height);
        assert(x >= 0);
        assert(y >= 0);
        write_color(ent.second, cold_val, hot_val, pixels + (y * width + x) * 3);
    }
    Magick::Image img(width, height, "RGB", MagickCore::StorageType::CharPixel, pixels);
    img.write(fn);
    delete [] pixels;
}

int main(int argc, char **argv) {
    Magick::InitializeMagick(*argv);
    if (argc < 9) {
        std::cout << "Usage: " << argv[0] << " <data file> <out file> <lat1> <lon1> <lat2> <lon2> <width> <height>\n";
        std::cout << "  lat1, lon1 and lat2, lon2 are the bounds of the output image\n";
        std::cout << "  height and width are the dimensions of the output image\n\n";
        std::cout << "Let t_i be the ith track, and p_i_j_lat, p_i_j_lon be the jth point of track i.\n";
        std::cout << "All coordinates are in floating point degrees\n";
        std::cout << "The data file is formatted as follows, where t_i is the ith track and\n";
        std::cout << "p_i_j is the jth point from track i:\n\n";
        std::cout << "<number of tracks (n)>\n";
        std::cout << "<number of points in t_1 (m)>\n";
        std::cout << "<p_1_1_lat p_1_1_lon>\n";
        std::cout << "...\n";
        std::cout << "<p_1_m_lat p_1_m_lon>\n";
        std::cout << "<number of points in t_2>\n";
        std::cout << "...\n";
        std::cout << "<number of points in t_n>\n";
        std::cout << "...\n";
        return 1;
    }
    std::vector<std::vector<Point>> tracks;
    int track_cnt;
    std::fstream input(argv[1]);
    input >> track_cnt;
    for (int i = 0; i < track_cnt; ++i) {
        int point_cnt;
        input >> point_cnt;
        std::vector<Point> track;
        track.reserve(point_cnt);
        for (int j = 0; j < point_cnt; ++j) {
            double lat_d, lon_d;
            input >> lat_d >> lon_d;
            track.push_back(Point(lat_d, lon_d, false));
        }
        tracks.push_back(std::move(track));
    }
    input.close();

    double lat1 = atof(argv[3]);
    double lon1 = atof(argv[4]);
    double lat2 = atof(argv[5]);
    double lon2 = atof(argv[6]);
    int width = atoi(argv[7]);
    int height = atoi(argv[8]);
    Point p1(lat1, lon1, false);
    Point p2(lat2, lon2, false);
    draw(argv[2], tracks, p1, p2, width, height);
}
