#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <initializer_list>
#include <limits>
#include <vector>
#include <tuple>
#include <memory>

using std::string;
using std::cout;
using std::cerr;
using std::endl;
using std::ofstream;
using std::vector;
using std::tuple;
using std::tie;
using std::make_tuple;
using std::shared_ptr;
using std::make_shared;
using std::dynamic_pointer_cast;

#include "io.h"
#include "matrix.h"

#include "MyObject.h"

Image
write_im(const vector<tuple<BinImage, uint, uint>> &comp, const Image &in)
{
	uint sx, sy;
	BinImage bin;
	Image dst(in.n_rows, in.n_cols);
	for(uint i = 0; i < dst.n_rows; ++i)
		for(uint j = 0; j < dst.n_cols; ++j)
			dst(i, j) = make_tuple(0, 0, 0);
	for(uint i = 0; i < comp.size(); ++i){
		tie(bin, sx, sy) = comp[i];
		for(uint k = 0; k < bin.n_rows; ++k)
			for(uint j = 0; j < bin.n_cols; ++j)
				if(bin(k, j) != 0){
					if(i == comp.size() - 1)
						dst(k+sx, j+sy) = make_tuple(255, 0, 0);
					else
						dst(k+sx, j+sy) = make_tuple(0, 255, 0);
				}
	}
	return dst.deep_copy();
}

BinImage binary(const Image &in)
{
	BinImage bin(in.n_rows, in.n_cols);
	uint i, j, r, g, b;
	for(i = 0; i < in.n_rows; ++i)
		for(j = 0; j < in.n_cols; ++j){
			tie(r, g, b) = in(i, j);
			if((r < 80) && (g < 80) && (b < 80))
				bin(i, j) = 0;
			else
				bin(i, j) = 1;
		}
	return bin.deep_copy();
}

void recursive_sv(BinImage &in, uint i, uint j)
{
	if(in(i, j) != 1)
		return;
	else{
		in(i, j) = 2;
		if(j < in.n_cols - 1)
			recursive_sv(in, i, j+1);
		if(j > 0)
			recursive_sv(in, i, j-1);
		if(i < in.n_rows - 1)
			recursive_sv(in, i+1, j);
		if(i > 0)
			recursive_sv(in, i-1, j);
	}
}

tuple<BinImage, uint, uint> comp_sv(BinImage &in, uint i, uint j)
{
	uint max_i = i, max_j = j, min_i = i, min_j = j;
	recursive_sv(in, i, j);
	bool begin_str = true;
	for(uint i1 = i; i1 < in.n_rows; ++i1){
		begin_str = true;
		for(uint j1 = 0; j1 < in.n_cols; ++j1){
			if(in(i1, j1) == 2){
				begin_str = false;
				if(i1 > max_i)
					max_i = i1;
				else if(i1 < min_i)
					min_i = i1;
				if(j1 > max_j)
					max_j = j1;
				else if(j1 < min_j)
					min_j = j1;
			}
		}
		if(begin_str)
			break;
	}
	BinImage subm =
		in.submatrix(min_i, min_j, max_i-min_i+1, max_j-min_j+1).deep_copy();
	for(uint i1 = 0; i1 < subm.n_rows; ++i1)
		for(uint j1 = 0; j1 < subm.n_cols; ++j1)
			if(subm(i1, j1) != 2)
				subm(i1, j1) = 0;
			else
				subm(i1, j1) = -1;
	for(uint i1 = 0; i1 < in.n_rows; ++i1)
		for(uint j1 = 0; j1 < in.n_cols; ++j1)
			if(in(i1, j1) == 2)
				in(i1, j1) = -1;
	return make_tuple(subm.deep_copy(), min_i, min_j);
}

vector<tuple<BinImage, uint, uint>> components(const BinImage &bin_image)
{
	auto comp_vec = vector<tuple<BinImage, uint, uint>>();
	BinImage comp_im = bin_image.deep_copy();
	for(uint i = 0; i < bin_image.n_rows; ++i)
		for(uint j = 0; j < bin_image.n_cols; ++j){
			if(comp_im(i, j) == 1){
				comp_vec.push_back(comp_sv(comp_im, i, j));
			}
		}
	return comp_vec;
}

bool on_gran(BinImage &in, uint i, uint j)
{
	return ((in(i,j) == -1) && ((i==0) || (j==0) || (i == in.n_rows - 1) ||
		(j == in.n_cols - 1) || (in(i-1,j) == 0) || (in(i+1,j) == 0) ||
		(in(i, j+1) == 0) || (in(i, j-1) == 0)));
}

tuple<uint, uint, uint, uint>
dist_t(BinImage &in)
{
	vector<tuple<uint, uint>> gran;
	uint min, max, min_r, max_r, rad, ci, cj, a, b;
	max = min_r = max_r = rad = ci = cj = min = a = b = 0;
	for(uint i = 0; i < in.n_rows; ++i)
		for(uint j = 0; j < in.n_cols; ++j)
			if(on_gran(in, i, j)){
				in(i, j) = -2;
				gran.push_back(make_tuple(i, j));
			}
	for(uint i = 0; i < in.n_rows; ++i)
		for(uint j = 0; j < in.n_cols; ++j)
			if(in(i, j) == -1){
				min = in.n_cols+in.n_rows;
				max = 0;
				for(uint k = 0; k < gran.size(); ++k){
					tie(a, b) = gran[k];
					rad = sqrt((a-i)*(a-i)+(b-j)*(b-j));
					if(rad < min)
						min = rad;
					if(rad > max)
						max = rad;
				}
				if(min > min_r){
					min_r = min;
					ci = i;
					cj = j;
					max_r = max;
				}
			}
	return make_tuple(ci, cj, min_r, max_r);
}

tuple<uint, uint, uint, uint>
find_center(tuple<BinImage, uint, uint> &vec)
{
	BinImage bin;
	uint si, sj, min_r, max_r, x, y;
	tie(bin, si, sj) = vec;
	tie(x, y, min_r, max_r) = dist_t(bin);
	return make_tuple(x+si, y+sj, min_r, max_r); 
}

uint square(const tuple<BinImage, uint, uint> &vec)
{
	uint s = 0, x, y;
	BinImage bin;
	tie(bin, x, y) = vec;
	for(uint i = 0; i < bin.n_rows; i++)
		for(uint j = 0; j < bin.n_cols; j++)
			if(bin(i,j) != 0)
				++s;
	return s;
}

uint perimetr(const tuple<BinImage, uint, uint> &vec)
{
	uint p = 0, x, y;
	BinImage bin;
	tie(bin, x, y) = vec;
	for(uint i = 0; i < bin.n_rows; i++)
		for(uint j = 0; j < bin.n_cols; j++)
			if(bin(i, j) == -2)
					++p;
	return p;
}

float compact(const tuple<BinImage, uint, uint> &vec)
{
	int s = square(vec), p = perimetr(vec);
	return (p*p)/s;
}

void add_axis(vector<shared_ptr<IObject>> &obj, uint x, uint y)
{
	obj.push_back(make_shared<Axis>(make_tuple(y,x)));
}

bool gear_br(const tuple<BinImage, uint, uint> &vec, uint x, uint y, uint nc)
{
	uint sx, sy, xc = 0, yc = 0, s = square(vec);
	BinImage bin;
	tie(bin, sx, sy) = vec;
	x -= sx;
	y -= sy;
	for(uint i = 0; i < bin.n_rows; ++i)
		for(uint j = 0; j < bin.n_cols; ++j)
			if(bin(i, j) == -1){
				xc += i;
				yc += j;
			}
	int X = xc / s + 4, Y = yc / s + 4;
	if (((X - x)*(X-x) <= 9) && ((Y-y)*(Y-y) <= 9) && (nc > 4))
		return false;
	else
		return true;
}

uint gear_nc
(tuple<BinImage, uint, uint> &bin, uint x, uint y, uint min_r, uint max_r)
{
	BinImage bin_n, bin_s;
	uint si = 0, sj = 0;
	tie(bin_n, si, sj) = bin;
	bin_s = bin_n.deep_copy();
	x -= si;
	y -= sj;
	for(uint i = 0; i < bin_s.n_rows; ++i)
		for(uint j = 0; j < bin_s.n_cols; ++j)
			if(sqrt((i-x)*(i-x)+(j-y)*(j-y)) < (min_r + max_r)/2)
				bin_s(i, j) = 0;
			else if (bin_s(i, j) != 0)
				bin_s(i, j) = 1;
	vector<tuple<BinImage, uint, uint>> comp_images = components(bin_s);
	return comp_images.size();
}

void add_gear(
	vector<shared_ptr<IObject>> &obj, uint x, uint y,
	tuple<BinImage, uint, uint> &bin, uint min_r, uint max_r)
{
	obj.push_back(make_shared<Gear>(make_tuple(y,x)));
	dynamic_pointer_cast<Gear>(obj[obj.size() - 1]) -> min_r = min_r+1;
	dynamic_pointer_cast<Gear>(obj[obj.size() - 1]) -> max_r = max_r;
	uint nc = gear_nc(bin, x, y, min_r, max_r); 
	dynamic_pointer_cast<Gear>(obj[obj.size() - 1]) -> num_cogs = nc;
	dynamic_pointer_cast<Gear>(obj[obj.size() - 1]) -> is_broken
		= gear_br(bin, x, y, nc);
}

vector<shared_ptr<IObject>>
classify(vector<tuple<BinImage, uint, uint>> &vec)
{
	vector<shared_ptr<IObject>> obj;
	uint x, y, s, min_r, max_r;
	for(uint i = 0; i < vec.size(); ++i){
		tie(x, y, min_r, max_r) = find_center(vec[i]);
		s = compact(vec[i]);
		if(s < 13)
			add_axis(obj, x, y);
		else
			add_gear(obj, x, y, vec[i], min_r, max_r);
	}
	return obj;
}

tuple<bool, uint, uint>
find_axis(vector<shared_ptr<IObject>> &vec)
{
	uint x, y;
	for(uint i = 0; i < vec.size(); ++i)
		if(dynamic_pointer_cast<Axis>(vec[i])){
			tie(x, y) = vec[i] -> location;
			return make_tuple(true, x, y);
		}
	return make_tuple(false, 0, 0);
}

tuple<bool, uint, uint>
find_broken(vector<shared_ptr<IObject>> &vec)
{
	uint x, y;
	for(uint i = 0; i < vec.size(); ++i)
		if(!dynamic_pointer_cast<Axis>(vec[i]) &&
		dynamic_pointer_cast<Gear>(vec[i]) -> is_broken){
			tie(x, y) = vec[i] -> location;
			return make_tuple(true, x, y);
		}
	return make_tuple(false, 0, 0);
}

bool
check(uint x, uint y, vector<shared_ptr<IObject>> &vec, uint min_r, uint max_r)
{
	uint x1, y1, rad, count = 0;
	for(uint i = 0; i < vec.size(); ++i){
		tie(x1, y1) = vec[i] -> location;
		if((x != x1) && (y != y1)){
			rad = sqrt((x-x1)*(x-x1)+(y-y1)*(y-y1));
			if((rad < dynamic_pointer_cast<Gear>(vec[i])->min_r + max_r) ||
				(rad < dynamic_pointer_cast<Gear>(vec[i])->max_r + min_r))
					return false;
			if(rad > dynamic_pointer_cast<Gear>(vec[i])->max_r + max_r)
				++count;
		}
	}
	if(count == vec.size() - 1)
		return false;
	else
		return true;
}

tuple<bool, BinImage, uint, uint>
what_image(vector<shared_ptr<IObject>> &vec, const char *name, uint n)
{
	char *str = new char[11];
	sprintf(str, "%c%c%c%c_%d.bmp", name[0], name[1], name[2], name[3], n);
	BinImage bin = binary(load_image(str));
	for(uint i = 0; i < bin.n_rows; ++i)
		for(uint j = 0; j < bin.n_cols; ++j)
			bin(i, j) = -bin(i,j);
	tuple<BinImage, uint, uint> v = make_tuple(bin, 0, 0);
	int x, y, min_r, max_r, xa = 0, ya = 0, xb = 0, yb = 0;
	bool ax = false, isb = false;
	tie(x, y, min_r, max_r) = find_center(v);
	tie(ax, xa, ya) = find_axis(vec);
	tie(isb, xb, yb) = find_broken(vec);
	if(ax)
		ax = check(xa, ya, vec, min_r, max_r);
	if(isb)
		isb = check(xb, yb, vec, min_r, max_r);
	if(ax)
		return make_tuple(true, bin.deep_copy(), xa-x, ya-y);
	else if (isb)
		return make_tuple(true, bin.deep_copy(), xb-x, yb-y);
	else
		return make_tuple(false, bin.deep_copy(), 0, 0);
}

tuple<uint, tuple<BinImage, uint, uint>>
	choose(vector<shared_ptr<IObject>> &vec, const char *name)
{
	bool b;
	BinImage bin;
	uint si, sj;
	for(uint i = 1; i <= 3; ++i){
		tie(b, bin, si, sj) = what_image(vec, name, i);
		if(b)
			return make_tuple(i, make_tuple(bin.deep_copy(), sj, si));
	}
	return make_tuple(0, make_tuple(bin.deep_copy(), 0, 0));
}

tuple<int, vector<shared_ptr<IObject>>, Image>
repair_mechanism(const Image& in, const char *name)
{
    // Base: return array of found objects and index of the correct gear
    // Bonus: return additional parameters of gears
    auto object_array = vector<shared_ptr<IObject>>();
    int result_idx = 0;
	BinImage bin = binary(in);
	tuple<BinImage, uint, uint>  finded;
	vector<tuple<BinImage, uint, uint>> comp_images = components(bin);
	object_array = classify(comp_images);
	tie(result_idx, finded) = choose(object_array, name);
	comp_images.push_back(finded);
	Image dst = write_im(comp_images, in);
    return make_tuple(result_idx, object_array, dst.deep_copy());
}

int main(int argc, char **argv)
{
    if (argc != 4)
    {
        cout << "Usage: " << endl << argv[0]
             << " <in_image.bmp> <out_image.bmp> <out_result.txt>" << endl;
        return 0;
    }

    try {
        Image src_image = load_image(argv[1]);
        ofstream fout(argv[3]);

        vector<shared_ptr<IObject>> object_array;
        Image dst_image;
        int result_idx;
        tie(result_idx, object_array, dst_image) =
			repair_mechanism(src_image, argv[1]);
        save_image(dst_image, argv[2]);

        fout << result_idx << endl;
        fout << object_array.size() << endl;
        for (const auto &obj : object_array)
            obj->Write(fout);

    } catch (const string &s) {
        cerr << "Error: " << s << endl;
        return 1;
    }
}
