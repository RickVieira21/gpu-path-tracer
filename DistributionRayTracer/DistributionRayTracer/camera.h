#ifndef CAMERA_H
#define CAMERA_H

#include <cmath>
#include <stdio.h>
using namespace std;

#include "vector.h"
#include "ray.h"
#include "maths.h"

class Camera
{

private:
	Vector eye, at, up; 
	float fovy, vnear, vfar, plane_dist, focal_ratio, aperture;
	float w, h;
	int res_x, res_y;
	Vector u, v, n;

public:

	Vector GetEye() { return eye; }
	int GetResX()  { return res_x; }
    int GetResY()  { return res_y; }
	float GetFov() { return fovy; }
	float GetPlaneDist() { return plane_dist; }
	float GetFar() {return vfar; }
	float GetAperture() { return aperture; }
	
    Camera( Vector from, Vector At, Vector Up, float angle, float hither, float yon, int ResX, int ResY, float Aperture_ratio, float Focal_ratio) {
	    eye = from;
	    at = At;
	    up = Up;
	    fovy = angle;
	    vnear = hither;
	    vfar = yon;
	    res_x = ResX;
	    res_y = ResY;
		focal_ratio = Focal_ratio;

        // set the camera frame uvn
        n = ( eye - at );
        plane_dist = n.length();
	    n = n / plane_dist;

	    u = up % n;
	    u = u / u.length();

	    v = n % u;

        //Dimensions of the vis window
	    h = 2 * plane_dist * tan( (PI * angle / 180) / 2.0f );
        w = ( (float) res_x / res_y ) * h;  

		aperture = Aperture_ratio * (w / res_x); //Lens aperture (diameter) = aperture_ratio * pixel_size (1 pixel=lente de raio 0.5)

		printf("\nwidth=%f height=%f fov=%f, viewplane distance=%f, pixel size=%.3f\n", w,h, fovy,plane_dist, w/res_x);
		if (Aperture_ratio != 0) printf("\nLens aperture = %.1f\n", Aperture_ratio);
    }

	void SetEye(Vector from) { 
		eye = from;
		// set the camera frame uvn
		n = (eye - at);
		plane_dist = n.length();
		n = n / plane_dist;
		u = up % n;
		u = u / u.length();
		v = n % u;
	}

	// Gera um raio primário que vai do olho (câmera) até um ponto de amostra no pixel na viewport
	Ray PrimaryRay(const Vector& pixel_sample) 
	{
		Vector ray_dir;  

		// Normalizar para colocar o pixel no sistema de coordenadas da câmera, no plano da imagem
		float px = (pixel_sample.x / res_x - 0.5f) * w; 
		float py = (pixel_sample.y / res_y - 0.5f) * h; 

		// Calcula a direção do raio no sistema de coordenadas da câmera:
		// 'u' e 'v' são os vetores base do plano da imagem (horizontal e vertical),
		// 'n' é o vetor normal para o plano da imagem (direção da câmera),
		// 'plane_dist' é a distância do plano da imagem ao olho
		ray_dir = (u * px) + (v * py) - (n * plane_dist);

		ray_dir.normalize();

		// Retorna o raio com origem no olho (posição da câmera) e direção calculada
		return Ray(eye, ray_dir);
	}

	// Gera um raio primário para câmeras do tipo "depth of field" 
	// O raio começa numa amostra na lente e vai até uma amostra no pixel na viewport
	Ray PrimaryRay(const Vector& lens_sample, const Vector& pixel_sample)
	{
		// Normaliza as coordenadas do pixel para o intervalo [-0.5, 0.5]
		float px = (pixel_sample.x / res_x - 0.5f) * w;
		float py = (pixel_sample.y / res_y - 0.5f) * h;

		// Calcula a posição do pixel no mundo, no plano da imagem, relativo à posição do olho
		Vector pixel_world = eye + (u * px) + (v * py) - (n * plane_dist);

		// Origem do raio é um ponto amostrado na lente da câmera (simulando abertura)
		// 'lens_sample' representa uma amostra aleatória dentro do disco da lente em coordenadas u e v
		Vector origin = eye + (u * lens_sample.x) + (v * lens_sample.y);

		// Direção do raio é do ponto da lente até o ponto do pixel no plano da imagem
		Vector dir = pixel_world - origin;

		dir.normalize();

		// Retorna o raio com origem no ponto amostrado na lente e direção calculada
		return Ray(origin, dir);
	}

};

#endif