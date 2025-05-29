/**
 * common.glsl
 * Common types and functions used for ray tracing.
 */

const float pi = 3.14159265358979;
const float epsilon = 0.001;

struct Ray {
    vec3 o;     // origin
    vec3 d;     // direction - always set with normalized vector
    float t;    // time, for motion blur
};

Ray createRay(vec3 o, vec3 d, float t)
{
    Ray r;
    r.o = o;
    r.d = d;
    r.t = t;
    return r;
}

Ray createRay(vec3 o, vec3 d)
{
    return createRay(o, d, 0.0);
}

vec3 pointOnRay(Ray r, float t)
{
    return r.o + r.d * t;
}

float gSeed = 0.0;

uint baseHash(uvec2 p)
{
    p = 1103515245U * ((p >> 1U) ^ (p.yx));
    uint h32 = 1103515245U * ((p.x) ^ (p.y>>3U));
    return h32 ^ (h32 >> 16);
}

float hash1(inout float seed) {
    uint n = baseHash(floatBitsToUint(vec2(seed += 0.1,seed += 0.1)));
    return float(n) / float(0xffffffffU);
}

vec2 hash2(inout float seed) {
    uint n = baseHash(floatBitsToUint(vec2(seed += 0.1,seed += 0.1)));
    uvec2 rz = uvec2(n, n * 48271U);
    return vec2(rz.xy & uvec2(0x7fffffffU)) / float(0x7fffffff);
}

vec3 hash3(inout float seed)
{
    uint n = baseHash(floatBitsToUint(vec2(seed += 0.1, seed += 0.1)));
    uvec3 rz = uvec3(n, n * 16807U, n * 48271U);
    return vec3(rz & uvec3(0x7fffffffU)) / float(0x7fffffff);
}

float rand(vec2 v)
{
    return fract(sin(dot(v.xy, vec2(12.9898, 78.233))) * 43758.5453);
}

vec3 toLinear(vec3 c)
{
    return pow(c, vec3(2.2));
}

vec3 toGamma(vec3 c)
{
    return pow(c, vec3(1.0 / 2.2));
}

vec2 randomInUnitDisk(inout float seed) {
    vec2 h = hash2(seed) * vec2(1.0, 6.28318530718);
    float phi = h.y;
    float r = sqrt(h.x);
	return r * vec2(sin(phi), cos(phi));
}

vec3 randomInUnitSphere(inout float seed)
{
    vec3 h = hash3(seed) * vec3(2.0, 6.28318530718, 1.0) - vec3(1.0, 0.0, 0.0);
    float phi = h.y;
    float r = pow(h.z, 1.0/3.0);
	return r * vec3(sqrt(1.0 - h.x * h.x) * vec2(sin(phi), cos(phi)), h.x);
}

vec3 randomUnitVector(inout float seed) //to be used in diffuse reflections with distribution cosine
{
    return(normalize(randomInUnitSphere(seed)));
}

struct Camera
{
    vec3 eye;
    vec3 u, v, n;
    float width, height;
    float lensRadius;
    float planeDist, focusDist;
    float time0, time1;
};

Camera createCamera(
    vec3 eye,
    vec3 at,
    vec3 worldUp,
    float fovy,
    float aspect,
    float aperture,  //diametro em multiplos do pixel size
    float focusDist,  //focal ratio
    float time0,
    float time1)
{
    Camera cam;
    if(aperture == 0.0) cam.focusDist = 1.0; //pinhole camera then focus in on vis plane
    else cam.focusDist = focusDist;
    vec3 w = eye - at;
    cam.planeDist = length(w);
    cam.height = 2.0 * cam.planeDist * tan(fovy * pi / 180.0 * 0.5);
    cam.width = aspect * cam.height;

    cam.lensRadius = aperture * 0.5 * cam.width / iResolution.x;  //aperture ratio * pixel size; (1 pixel=lente raio 0.5)
    cam.eye = eye;
    cam.n = normalize(w);
    cam.u = normalize(cross(worldUp, cam.n));
    cam.v = cross(cam.n, cam.u);
    cam.time0 = time0;
    cam.time1 = time1;
    return cam;
}

// done
Ray getRay(Camera cam, vec2 pixel_sample)  //rnd pixel_sample viewport coordinates
{
    vec2 ls = cam.lensRadius * randomInUnitDisk(gSeed);  //ls - lens sample for DOF
    float time = cam.time0 + hash1(gSeed) * (cam.time1 - cam.time0);
    
    //Calculate eye_offset and ray direction
    vec3 offset = cam.u * ls.x + cam.v * ls.y;
    vec3 eye_offset = cam.eye + offset;

    // Calcula o ponto alvo no plano de foco usando as coordenadas do pixel_sample
    vec3 target = cam.eye - cam.n * cam.focusDist + cam.u * (pixel_sample.x * cam.width)
                + cam.v * (pixel_sample.y * cam.height);

    vec3 ray_direction = normalize(target - eye_offset);

    return createRay(eye_offset, normalize(ray_direction), time);
}

// MT_ material type
#define MT_DIFFUSE 0
#define MT_METAL 1
#define MT_DIELECTRIC 2

struct Material
{
    int type;
    vec3 albedo;  //diffuse color
    vec3 specColor;  //the color tint for specular reflections. for metals and opaque dieletrics like coloured glossy plastic
    vec3 emissive; //
    float roughness; // controls roughness for metals. It can be used for rough refractions
    float refIdx; // index of refraction for Dielectric
    vec3 refractColor; // absorption for beer's law
};

Material createDiffuseMaterial(vec3 albedo)
{
    Material m;
    m.type = MT_DIFFUSE;
    m.albedo = albedo;
    m.specColor = vec3(0.0);
    m.roughness = 1.0;  //ser usado na iluminação direta
    m.refIdx = 1.0;
    m.refractColor = vec3(0.0);
    m.emissive = vec3(0.0);
    return m;
}

Material createMetalMaterial(vec3 specClr, float roughness)
{
    Material m;
    m.type = MT_METAL;
    m.albedo = vec3(0.0);
    m.specColor = specClr;
    m.roughness = roughness;
    m.emissive = vec3(0.0);
    return m;
}

Material createDielectricMaterial(vec3 refractClr, float refIdx, float roughness)
{
    Material m;
    m.type = MT_DIELECTRIC;
    m.albedo = vec3(0.0);
    m.specColor = vec3(0.04);
    m.refIdx = refIdx;
    m.refractColor = refractClr;  
    m.roughness = roughness;
    m.emissive = vec3(0.0);
    return m;
}

struct HitRecord
{
    vec3 pos;
    vec3 normal;
    float t;            // ray parameter
    Material material;
};

// done
float schlick(float cosine, float refIdx)
{
    float ior1 = 1.0;
    float reflect0 = (ior1 - refIdx) / (ior1 + refIdx);
    float r0 = reflect0 * reflect0;
    return r0 + (1.0 - r0) * pow(1.0 - cosine, 5.0);
}


bool scatter(Ray rIn, HitRecord rec, out vec3 atten, out Ray rScattered)
{
    if (rec.material.type == MT_DIFFUSE)
    {
        // Lambertian scatter
        vec3 scatterDirection = rec.normal + randomUnitVector(gSeed);

        // Corrige caso o vetor seja quase zero
        if (length(scatterDirection) < 1e-8)
            scatterDirection = rec.normal;

        rScattered = createRay(rec.pos, scatterDirection);  // <-- usa createRay
        atten = rec.material.albedo * max(dot(normalize(rScattered.d), rec.normal), 0.0) / 3.141592;
        return true;
    }

    if (rec.material.type == MT_METAL)
    {
        vec3 reflected = reflect(normalize(rIn.d), rec.normal);
        vec3 scatteredDir = reflected + rec.material.roughness * randomInUnitSphere(gSeed);

        rScattered = createRay(rec.pos, scatteredDir);
        atten = rec.material.specColor;

        return dot(rScattered.d, rec.normal) > 0.0;
    }

    if (rec.material.type == MT_DIELECTRIC)
    {
        atten = vec3(1.0);
        vec3 outwardNormal;
        float niOverNt;
        float cosine;

        vec3 rdNorm = normalize(rIn.d);

        if (dot(rdNorm, rec.normal) > 0.0)
        {
            outwardNormal = -rec.normal;
            niOverNt = rec.material.refIdx;
            cosine = rec.material.refIdx * dot(rdNorm, rec.normal);

            float dist = length(rIn.d); // distância percorrida (opcional, se quiseres ajustar)
            atten = exp(-rec.material.refractColor * dist);
        }
        else
        {
            outwardNormal = rec.normal;
            niOverNt = 1.0 / rec.material.refIdx;
            cosine = -dot(rdNorm, rec.normal);
        }

        vec3 refracted = refract(rdNorm, outwardNormal, niOverNt);
        bool canRefract = length(refracted) > 0.0001;

        float reflectProb = canRefract ? schlick(cosine, rec.material.refIdx) : 1.0;

        if (hash1(gSeed) < reflectProb)
        {
            vec3 reflected = reflect(rdNorm, rec.normal);
            rScattered = createRay(rec.pos, reflected);
        }
        else
        {
            rScattered = createRay(rec.pos, refracted);
        }

        return true;
    }

    return false;
}


struct Triangle {vec3 a; vec3 b; vec3 c; };

Triangle createTriangle(vec3 v0, vec3 v1, vec3 v2)
{
    Triangle t;
    t.a = v0; t.b = v1; t.c = v2;
    return t;
}

// done
bool hit_triangle(Triangle t, Ray r, float tmin, float tmax, out HitRecord rec)
{
    vec3 v0 = t.a;
    vec3 v1 = t.b;
    vec3 v2 = t.c;

    vec3 edge1 = v1 - v0;
    vec3 edge2 = v2 - v0;

    vec3 h = cross(r.d, edge2);
    float a = dot(edge1, h);

    //if (abs(a) < 1e-8)
    //    return false; // Raio paralelo ao triângulo

    float f = 1.0 / a;
    vec3 s = r.o - v0;
    float u = f * dot(s, h);

    if (u < 0.0 || u > 1.0)
        return false;

    vec3 q = cross(s, edge1);
    float v = f * dot(r.d, q);

    if (v < 0.0 || u + v > 1.0)
        return false;

    float t_intersect = f * dot(edge2, q);

    if (t_intersect < tmin || t_intersect > tmax)
        return false;

    // Interseção válida
    rec.t = t_intersect;
    rec.pos = r.o + t_intersect * r.d;
    rec.normal = normalize(cross(edge1, edge2));

    return true;
}



struct Quad {vec3 a; vec3 b; vec3 c; vec3 d; };

Quad createQuad(vec3 v0, vec3 v1, vec3 v2, vec3 v3)
{
    Quad q;
    q.a = v0; q.b = v1; q.c = v2; q.d = v3;
    return q;
}

bool hit_quad(Quad q, Ray r, float tmin, float tmax, out HitRecord rec)
{
    if(hit_triangle(createTriangle(q.a, q.b, q.c), r, tmin, rec.t, rec)) return true;
    else if(hit_triangle(createTriangle(q.a, q.c, q.d), r, tmin, rec.t, rec)) return true;
    else return false;  
}


struct Sphere
{
    vec3 center;
    float radius;
};

Sphere createSphere(vec3 center, float radius)
{
    Sphere s;
    s.center = center;
    s.radius = radius;
    return s;
}


struct MovingSphere
{
    vec3 center0, center1;
    float radius;
    float time0, time1;
};

MovingSphere createMovingSphere(vec3 center0, vec3 center1, float radius, float time0, float time1)
{
    MovingSphere s;
    s.center0 = center0;
    s.center1 = center1;
    s.radius = radius;
    s.time0 = time0;
    s.time1 = time1;
    return s;
}

// done
vec3 center(MovingSphere mvsphere, float time)
{
    float t = (time - mvsphere.time0) / (mvsphere.time1 - mvsphere.time0);
    return mix(mvsphere.center0, mvsphere.center1, clamp(t, 0.0, 1.0));
}


/*
 * The function naming convention changes with these functions to show that they implement a sort of interface for
 * the book's notion of "hittable". E.g. hit_<type>.
 */
// done
bool hit_sphere(Sphere s, Ray r, float tmin, float tmax, out HitRecord rec)
{
    vec3 oc = r.o - s.center;
    float a = dot(r.d, r.d);
    float b = dot(oc, r.d);
    float c = dot(oc, oc) - s.radius * s.radius;
    float discriminant = b * b - a * c;

    if (discriminant > 0.0) {
        float sqrtD = sqrt(discriminant);

        float temp = (-b - sqrtD) / a;
        if (temp < tmax && temp > tmin) {
            rec.t = temp;
            rec.pos = pointOnRay(r, rec.t);
            rec.normal = (rec.pos - s.center) / s.radius;

            // Face normal logic
            bool frontFace = dot(r.d, rec.normal) < 0.0;
            rec.normal = frontFace ? rec.normal : -rec.normal;

            return true;
        }
        temp = (-b + sqrtD) / a;

        if (temp < tmax && temp > tmin) {
            rec.t = temp;
            rec.pos = pointOnRay(r, rec.t);
            rec.normal = (rec.pos - s.center) / s.radius;

            // Face normal logic
            bool frontFace = dot(r.d, rec.normal) < 0.0;
            rec.normal = frontFace ? rec.normal : -rec.normal;

            return true;
        }
    }
    return false;
}

// done
bool hit_movingSphere(MovingSphere s, Ray r, float tmin, float tmax, out HitRecord rec)
{
    // 1. Interpolate the center at the ray's time
    float timeFraction = (r.t - s.time0) / (s.time1 - s.time0);
    vec3 center = mix(s.center0, s.center1, clamp(timeFraction, 0.0, 1.0));  // Linear interpolation

    // 2. Ray-sphere intersection (same as static sphere)
    vec3 oc = r.o - center;
    float a = dot(r.d, r.d);
    float b = dot(oc, r.d);
    float c = dot(oc, oc) - s.radius * s.radius;
    float discriminant = b * b - a * c;

    if (discriminant > 0.0) {
        float sqrtD = sqrt(discriminant);

        float temp = (-b - sqrtD) / a;
        if (temp < tmax && temp > tmin) {
            rec.t = temp;
            rec.pos = pointOnRay(r, rec.t);
            rec.normal = (rec.pos - center) / s.radius;

            // Flip normal if hitting from inside
            bool frontFace = dot(r.d, rec.normal) < 0.0;
            rec.normal = frontFace ? rec.normal : -rec.normal;

            return true;
        }

        temp = (-b + sqrtD) / a;
        if (temp < tmax && temp > tmin) {
            rec.t = temp;
            rec.pos = pointOnRay(r, rec.t);
            rec.normal = (rec.pos - center) / s.radius;

            bool frontFace = dot(r.d, rec.normal) < 0.0;
            rec.normal = frontFace ? rec.normal : -rec.normal;

            return true;
        }
    }
    return false;
}

struct pointLight {
    vec3 pos;
    vec3 color;
};

pointLight createPointLight(vec3 pos, vec3 color) 
{
    pointLight l;
    l.pos = pos;
    l.color = color;
    return l;
}