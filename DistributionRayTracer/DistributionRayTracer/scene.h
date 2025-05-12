#ifndef SCENE_H
#define SCENE_H

#include <vector>
#include <cmath>
#include <IL/il.h>
using namespace std;

#include "camera.h"
#include "color.h"
#include "vector.h"
#include "ray.h"
#include "boundingBox.h"

//Light types
typedef enum {PUNCTUAL, QUAD} lightType;

//Skybox images constant symbolics
typedef enum { RIGHT, LEFT, TOP, BOTTOM, FRONT, BACK } CubeMap;

//Type of acceleration structure
typedef enum { NONE, GRID_ACC, BVH_ACC }  accelerator;

struct HitRecord
{
	bool isHit = false;
	Vector normal;
	float t = FLT_MAX;            // ray parameter
	Vector bary;  //barycentric coordinates for interpolation in a triangle
	Vector texUV;  //interpolated texCoord in the hit point
};


class Material
{
public:
	
	Material() :
		m_diffColor(Color(0.2f, 0.2f, 0.2f)), m_Diff(0.2f), m_specColor(Color(1.0f, 1.0f, 1.0f)), m_Spec(0.8f), m_Shine(20), m_Refl(1.0f), m_T(0.0f), m_RIndex(1.0f) {};

	Material(Color& c, float Kd, Color& cs, float Ks, float Shine, float T, float ior) {
		m_diffColor = c; m_Diff = Kd; m_specColor = cs; m_Spec = Ks; m_Shine = Shine; m_Refl = Ks; m_T = T; m_RIndex = ior;
	}

	void SetDiffColor(Color& a_Color) { m_diffColor = a_Color; }
	Color GetDiffColor() { return m_diffColor; }
	void SetSpecColor(Color& a_Color) { m_specColor = a_Color; }
	Color GetSpecColor() { return m_specColor; }
	void SetDiffuse(float a_Diff) { m_Diff = a_Diff; }
	void SetSpecular(float a_Spec) { m_Spec = a_Spec; }
	void SetShine(float a_Shine) { m_Shine = a_Shine; }
	void SetReflection(float a_Refl) { m_Refl = a_Refl; }
	void SetTransmittance(float a_T) { m_T = a_T; }
	float GetSpecular() { return m_Spec; }
	float GetDiffuse() { return m_Diff; }
	float GetShine() { return m_Shine; }
	float GetReflection() { return m_Refl; }
	float GetTransmittance() { return m_T; }
	void SetRefrIndex(float a_ior) { m_RIndex = a_ior; }
	float GetRefrIndex() { return m_RIndex; }
private:
	Color m_diffColor, m_specColor;
	float m_Refl, m_T;
	float m_Diff, m_Shine, m_Spec;
	float m_RIndex;
};

class Light
{

private:
	Vector e1, e2; //Light frame in World coordinates
	

public:

	Vector position;
	Color emission;
	float area;
	Vector normal;
	lightType type;
	unsigned int gridRes;  // resolution of a regular grid; to be used ONLY without Antialiasing; otherwise it should be 0

	Light(Vector& pos, Color& col, Vector& v1, Vector& v2, unsigned int grid_res) {
		type = QUAD;
		position = pos;   //position of point light or the center in the area light
		emission = col;
		gridRes = grid_res;

		e1 = v1 - position;
		e2 = v2 - position;

		area = (e1 % e2).length();
		normal = (e1 % e2).normalize();
	}

	Light(Vector& pos, Color& col) {
		type = PUNCTUAL;
		position = pos;   //position of point light or the center in the area light
		emission = col;
	}

	Vector getAreaLightPoint(const Vector& sample)   //get a point in WC
	{
		return(position + e1 * sample.x + e2 * sample.y);
	}
};

class Object
{
public:

	Material* GetMaterial() { return m_Material; }
	void SetMaterial( Material *a_Mat ) { m_Material = a_Mat; }
	virtual HitRecord hit( Ray& r) = 0;
	virtual AABB GetBoundingBox() { return AABB(); }
	Vector getCentroid(void) { return GetBoundingBox().centroid(); }

protected:
	Material* m_Material;
	
};

class Plane : public Object
{
protected:
  Vector	 PN;
  float 	 D;

public:
		 Plane		(Vector& PNc, float Dc);
		 Plane		(Vector& P0, Vector& P1, Vector& P2);

		 HitRecord hit(Ray& r);
};

class Triangle : public Object
{
public:
	Triangle	(Vector& P0, Vector& P1, Vector& P2);
	AABB GetBoundingBox(void);
	HitRecord hit(Ray& r);
	
protected:
	Vector points[3];
	Vector n1, n2, n3;   //vertex normals from obj files
	Vector t1, t2, t3;   //vertex textures from obj files
	Vector Min, Max;
};



class Sphere : public Object
{
public:
	Sphere( Vector& a_center, float a_radius ) : center( a_center ), SqRadius( a_radius * a_radius ), radius( a_radius ) {};
	HitRecord hit(Ray& r);
	AABB GetBoundingBox(void);

private:
	Vector center;
	float radius, SqRadius;
};

class aaBox : public Object   //Axis aligned box: another geometric object
{
public:
	aaBox(Vector& minPoint, Vector& maxPoint);
	AABB GetBoundingBox(void);
	HitRecord hit(Ray& r);

private:
	Vector min;
	Vector max;

	Vector Normal;
};


class Scene
{
public:
	Scene();
	virtual ~Scene();
	
	Camera* GetCamera() { return camera; }
	Color GetBackgroundColor() { return bgColor; }
	bool GetSkyBoxFlg() { return SkyBoxFlg; }
	Color GetSkyboxColor(Ray& r);
	unsigned int GetSamplesPerPixel() { return samples_per_pixel; }
	accelerator GetAccelStruct() { return accel_struc_type; }

	void SetBackgroundColor(Color a_bgColor) { bgColor = a_bgColor; }
	void SetSkyBoxFlg(bool a_skybox_flg) {SkyBoxFlg = a_skybox_flg;}
	void LoadSkybox(const char*);
	void SetCamera(Camera *a_camera) {camera = a_camera; }
	void SetAccelStruct(accelerator accel_t) { accel_struc_type = accel_t; }
	void SetSamplesPerPixel(unsigned int spp) { samples_per_pixel = spp; }
	
	int getNumObjects( );
	void addObject( Object* o );
	Object* getObject( unsigned int index );
	
	int getNumLights( );
	void addLight( Light* l );
	Light* getLight( unsigned int index );

	bool load_p3f(const char *name);  //Load NFF file method
	void create_random_scene();
	
private:
	vector<Object *> objects;
	vector<Light *> lights;

	Camera* camera;
	Color bgColor;  //Background color
	unsigned int samples_per_pixel;  // samples per pixel
	accelerator accel_struc_type; 

	bool SkyBoxFlg;
	struct {
		ILubyte *img;
		unsigned int resX;
		unsigned int resY;
		unsigned int BPP; //bytes per pixel
	} skybox_img[6];

};

#endif