#define _CRT_SECURE_NO_WARNINGS 1
#include <vector>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

#include <random>
#include <algorithm>

static std::default_random_engine engine(10);
static std::normal_distribution<double> normal(0, 1);

std::array<double, 3> random_dir() {
	double xDir = normal(engine);
	double yDir = normal(engine);
	double zDir = normal(engine);
	double norm = sqrt(xDir * xDir + yDir * yDir + zDir * zDir);
	xDir /= norm;
	yDir /= norm;
	zDir /= norm;
	return {xDir, yDir, zDir};
}

void sliced_matching(const std::vector<double>& image_source, const std::vector<double>& color_source, int W, int H, std::vector<double>& image_result) {
	// This is a copy 
	image_result = image_source;

	for(int iter = 0; iter < 100; ++iter) {
		auto [xDir, yDir, zDir] = random_dir();
		std::vector<std::pair<double, int>> result_sorted(W * H);
		std::vector<double> sortedTarget(W * H);
		for(int i = 0; i < W * H; ++i) {
			double dotResult = image_result[i * 3] * xDir + image_result[i * 3 + 1] * yDir + image_result[i * 3 + 2] * zDir;
			double dotTarget = color_source[i * 3] * xDir + color_source[i * 3 + 1] * yDir + color_source[i * 3 + 2] * zDir;
			result_sorted[i] = std::make_pair(dotResult, i);
			sortedTarget[i] = dotTarget;
		}
		std::sort(result_sorted.begin(), result_sorted.end());
		std::sort(sortedTarget.begin(), sortedTarget.end());

		for(int i = 0; i < W * H; ++i) {
			double motionAmount = sortedTarget[i] - result_sorted[i].first;
			int index = result_sorted[i].second;
			image_result[index * 3 + 0] += motionAmount * xDir;
			image_result[index * 3 + 1] += motionAmount * yDir;
			image_result[index * 3 + 2] += motionAmount * zDir;
		}
	}
}

unsigned char clamp(double x, double min, double max) {
	if(x < min) return min;
	if(x > max) return max;
	return x;
}

int main() {
	int W, H, C;
	
	// stbi_set_flip_vertically_on_load(false);
	unsigned char *image_source_ptr = stbi_load("8733654151_b9422bb2ec_k.jpg",
                                 &W,
                                 &H,
                                 &C,
                                 STBI_rgb);
	unsigned char *color_source_ptr = stbi_load("redim.jpg",
								 &W,
								 &H,
								 &C,
								 STBI_rgb);

	size_t total_size = W*H*3;
	std::vector<double> image_source(total_size);
	std::vector<double> color_source(total_size);
	for (int i = 0; i < total_size; ++i) {
		image_source[i] = (double)image_source_ptr[i];
		color_source[i] = (double)color_source_ptr[i];
	}

	std::vector<double> result(total_size);
	sliced_matching(image_source, color_source, W, H, result);
	
	unsigned char *image_result = new unsigned char[W*H*3];
	for (int i = 0; i < H; ++i) {
		for (int j = 0; j < W; ++j) {
			image_result[(i*W + j) * 3 + 0] = clamp(result[(i*W + j) * 3 + 0], 0, 255);
			image_result[(i*W + j) * 3 + 1] = clamp(result[(i*W + j) * 3 + 1], 0, 255);
			image_result[(i*W + j) * 3 + 2] = clamp(result[(i*W + j) * 3 + 2], 0, 255);
		}
	}
	stbi_write_png("image.png", W, H, 3, image_result, 0);

	return 0;
}