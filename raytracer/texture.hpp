#pragma once 

#include <vector>
#include "vector.h"

class Texture {
private:
    std::vector<Vector3> pixels;
    int width, height, channels;
public:
    Texture(unsigned char *pixels_src, int width, int height, int channels): width(width), height(height), channels(channels) {
        this->pixels = std::vector<Vector3>(width * height);
        for (int i = 0; i < width * height; ++i) {
            // NOTE: The initial picture is gamma corrected. We need to correct it back and also divide by 255
            this->pixels[i] = Vector3(pixels_src[i * 3 + 0], pixels_src[i * 3 + 1], pixels_src[i * 3 + 2]);
            this->pixels[i] = this->pixels[i] / 255.0;
            this->pixels[i] = reverse_gamma_correct(this->pixels[i]);
        }
    }

    Texture(const Vector3& color) {
        this->pixels = std::vector<Vector3>(1, color);
        this->width = 1;
        this->height = 1;
        this->channels = 3;
    }

    Vector3 get_pixel(Vector3 uv) const {
        int i = std::min(std::max((int)(uv.x * width ), 0), width - 1);
        int j = height - 1 - std::min(std::max((int)(uv.y * height), 0), height - 1);
        std::swap(i, j);
        return pixels[i * width + j];
    }
};