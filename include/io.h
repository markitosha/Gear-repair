#pragma once

#include "matrix.h"
#include "EasyBMP.h"

#include <tuple>

typedef Matrix<std::tuple<uint, uint, uint>> Image;

typedef Matrix<int> BinImage;

Image load_image(const char*);
void save_image(const Image&, const char*);
