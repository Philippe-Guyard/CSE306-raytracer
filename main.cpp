#include <vector>
#include <string>
#include <iostream>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"
 
#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

#include "raytracer/ray.h"
#include "raytracer/sphere.h"
#include "raytracer/vector.h"
#include "raytracer/scene.h"
#include "raytracer/triangles.h"

class Image {
private:
    int W, H;
    std::vector<unsigned char> data;
public:
    Image(int W, int H): W(W), H(H), data(W*H*3, 0) {};

    ~Image() {
        data.clear();
    };

    void set_pixel(int i, int j, const Vector3& color, bool gamma_correction = true) {
        Vector3 c = gamma_correction ? gamma_correct(color) : color;
        data[(i*W + j) * 3 + 0] = std::min(c.x, 255.0);
        data[(i*W + j) * 3 + 1] = std::min(c.y, 255.0);
        data[(i*W + j) * 3 + 2] = std::min(c.z, 255.0);
    }

    void save(const std::string& filename) {    
        stbi_write_png(filename.c_str(), W, H, 3, &data[0], 0);
    }
};

void box_muller(double &x, double &y) {
    double r1 = Rng::random_real();
    double r2 = Rng::random_real();
    x = sqrt(-2 * log(r1)) * cos(2 * M_PI * r2);
    y = sqrt(-2 * log(r1)) * sin(2 * M_PI * r2);
}

int main() {
    int W = 512;
    int H = 512;
    double alpha = 60.*M_PI/180.;

    LightSource light(Vector3(-10, 20, 40), 2E10);

    // Sphere center(Vector3(-10, -5, 0), 10, Vector3(0., 0.5, 1.), true, false);
    // Sphere center2(Vector3(10, -5, 0), 10, Vector3(1., 1., 1.), false, false);
    Sphere left_wall(Vector3(-1000, 0, 0), 940, Vector3(0.5, 0.8, 0.1));
	Sphere right_wall(Vector3(1000, 0, 0), 940, Vector3(0.9, 0.2, 0.3));
	Sphere ceiling(Vector3(0, 1000, 0), 940, Vector3(0.3, 0.5, 0.3));
	Sphere floor(Vector3(0, -955, 0), 940, Vector3(0.6, 0.5, 0.7));
	Sphere front_wall(Vector3(0, 0, -1000), 940, Vector3(0.1, 0.6, 0.7));
	Sphere behind_wall(Vector3(0, 0, 1000), 940, Vector3(0.0, 0.2, 0.9));
    Triangle triangle(Vector3(-5, 12, 0), Vector3(5, 12, 0), Vector3(0, 22, 0), Vector3(1, 1, 1));
    // Sphere random(Vector3(0, 6, 0), 6, Vector3(1, 1, 1));

    Scene scene(light);
    // scene.add_object(&center);
    // scene.add_object(&center2);
    scene.add_object(&left_wall);
    scene.add_object(&right_wall);
    scene.add_object(&ceiling);
    scene.add_object(&floor);
    scene.add_object(&front_wall);
    scene.add_object(&behind_wall);
    scene.add_object(&triangle);
    // scene.add_object(&random);

    Vector3 camera_center(0, 0, 55);
	
    Image img(W, H);
#pragma omp parallel for schedule(dynamic, 1)
    for (int i = 0; i < H; i++) {
        for (int j = 0; j < W; j++) {
            Vector3 color = Vector3(0, 0, 0);
            for(int k = 0; k < 32; k++) {
                double randomX, randomY;
                box_muller(randomX, randomY);
                Vector3 dir = Vector3(j - W / 2. + 0.5 + randomX * 0.5, 
                                     -i + H / 2. + 0.5 + randomY * 0.5, 
                                     -W/(2.*tan(alpha/2.)));
                Ray ray(camera_center, dir);
                color += scene.get_color(ray);
            }
            color = color / 32.;
            img.set_pixel(i, j, color);
        }
    }
    img.save("out.png");

    return 0;
}