#include <iostream>
#include <sstream>
#include <vector>
#include <memory>
#include <stdint.h>
#include <cmath>
#include <random>
#include <omp.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include <string>
#include <fstream>
#include <sstream>
#include <compare>
#include <filesystem>
#include <chrono>
#include <cstdlib>

#include <opencv2/opencv.hpp>

#define THREADS 10


using vec = std::vector<double>;
using matrix = std::vector<std::vector<double>>;
using mat3D = std::vector<std::vector<std::vector<double>>>;
using u32 = uint_fast32_t;
using i32 = int_fast32_t;
