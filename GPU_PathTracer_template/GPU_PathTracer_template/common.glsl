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

// Gera um raio da câmara para um ponto no plano de imagem, com profundidade de campo
Ray getRay(Camera cam, vec2 pixelSample) {
    // Sample aleatório na lente para DOF (profundidade de campo)
    vec2 lensSample = cam.lensRadius * randomInUnitDisk(gSeed);

    // Tempo aleatório para motion blur
    float rayTime = cam.time0 + hash1(gSeed) * (cam.time1 - cam.time0);

    // Distâncias relevantes
    float d = cam.planeDist;                // distância ao plano de imagem
    float f = d * cam.focusDist;            // distância ao plano de foco

    // Conversão das coordenadas do pixel para o plano de imagem 
    float px = (pixelSample.x + 0.5) / iResolution.x - 0.5;
    float py = (pixelSample.y + 0.5) / iResolution.y - 0.5;

    // Ponto no plano de imagem (em coordenadas da câmara)
    vec3 imagePlanePoint = vec3(px * cam.width, py * cam.height, -d);

    // Converte o ponto para o plano de foco 
    vec3 focalPoint = imagePlanePoint * cam.focusDist;

    // Deslocamento do olho (camera origin) devido à abertura da lente
    vec3 eyeOffset = cam.eye + lensSample.x * cam.u + lensSample.y * cam.v;

    // Direção do raio (do olho deslocado para o ponto de foco)
    vec3 rayDir = (focalPoint.x - lensSample.x) * cam.u
                + (focalPoint.y - lensSample.y) * cam.v
                - f * cam.n;

    return Ray(eyeOffset, normalize(rayDir), rayTime);
}


// MT_ material type
#define MT_DIFFUSE 0
#define MT_METAL 1
#define MT_DIELECTRIC 2
#define MT_PLASTIC 3  


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


Material createPlasticMaterial(vec3 albedo, float roughness)
{
    Material m;
    m.type = MT_PLASTIC;
    m.albedo = albedo;
    m.specColor = vec3(0.04); 
    m.roughness = roughness;
    m.refIdx = 1.0;
    m.refractColor = vec3(0.0);
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


//MICROFACETS

vec3 fresnelSchlick(float cosTheta, vec3 F0) {
    return F0 + (1.0 - F0) * pow(1.0 - cosTheta, 5.0);
}

float D_GGX(float NoH, float roughness) {
    float alpha = roughness * roughness;
    float alpha2 = alpha * alpha;
    float NoH2 = NoH * NoH;
    float b = (NoH2 * (alpha2 - 1.0) + 1.0);
    return alpha2 / (pi * b * b);
}

float G1_GGX_Schlick(float NoV, float roughness) {
    // float r = roughness; // original
    float r = 0.5 + 0.5 * roughness; // Disney remapping
    float k = (r * r) / 2.0;
    float denom = NoV * (1.0 - k) + k;
    return max(NoV, 0.001) / denom;
}

float G_Smith(float NoV, float NoL, float roughness) {
    float g1_l = G1_GGX_Schlick(NoL, roughness);
    float g1_v = G1_GGX_Schlick(NoV, roughness);
    return g1_l * g1_v;
}

vec3 sampleGGX(vec3 V, float roughness, inout float seed)
{
    float a = roughness * roughness;

    float phi = 2.0 * pi * hash1(seed);
    float cosTheta = sqrt((1.0 - hash1(seed)) / (1.0 + (a*a - 1.0) * hash1(seed)));
    float sinTheta = sqrt(1.0 - cosTheta * cosTheta);

    vec3 H = vec3(cos(phi) * sinTheta, sin(phi) * sinTheta, cosTheta);

    // Transformar H para estar alinhado com a direção V
    vec3 up = abs(V.z) < 0.999 ? vec3(0.0, 0.0, 1.0) : vec3(1.0, 0.0, 0.0);
    vec3 tangentX = normalize(cross(up, V));
    vec3 tangentY = cross(V, tangentX);

    return normalize(H.x * tangentX + H.y * tangentY + H.z * V);
}




bool scatter(Ray rIn, HitRecord rec, out vec3 atten, out Ray rScattered)
{
    if (rec.material.type == MT_DIFFUSE)
    {
        // Lambertian scatter
        vec3 scatterDirection = rec.normal + randomUnitVector(gSeed);

        if (length(scatterDirection) < 1e-8)
            scatterDirection = rec.normal;

        rScattered = createRay(rec.pos, scatterDirection);
        atten = rec.material.albedo * max(dot(normalize(rScattered.d), rec.normal), 0.0) / 3.141592;
        return true;
    }

    if (rec.material.type == MT_METAL)
    {
        // Microfacet reflection using GGX
        vec3 V = normalize(-rIn.d);      // View direction
        vec3 N = rec.normal;             // Normal
        float roughness = rec.material.roughness;

        // Sample half-vector H with GGX
        vec3 H = sampleGGX(N, roughness, gSeed);
        vec3 L = reflect(-V, H);         // Reflected direction

        // Valid reflection?
        if (dot(L, N) <= 0.0)
            return false;

        rScattered = createRay(rec.pos + N * epsilon, normalize(L));
        
        float NoH = max(dot(N, H), 0.001);
        float NoV = max(dot(N, V), 0.001);
        float NoL = max(dot(N, L), 0.001);

        vec3 F = fresnelSchlick(max(dot(H, V), 0.0), rec.material.specColor);
        float D = D_GGX(NoH, roughness);
        float G = G_Smith(NoV, NoL, roughness);

        vec3 numerator = F * D * G;
        float denominator = 4.0 * max(dot(N, V), 0.001) * max(dot(N, L), 0.001);
        atten = numerator / max(denominator, 0.001);

        return true;
    }

    if (rec.material.type == MT_PLASTIC)
    {
        vec3 V = normalize(-rIn.d);
        vec3 N = rec.normal;
        float roughness = rec.material.roughness;

        // Calcula Fresnel para ponderar a escolha probabilística
        float cosTheta = max(dot(N, V), 0.0);
        float fresnel = fresnelSchlick(cosTheta, rec.material.specColor).r;

        if (hash1(gSeed) < fresnel)  // Especular via microfacet GGX
        {
            vec3 H = sampleGGX(N, roughness, gSeed);
            vec3 L = reflect(-V, H);

            if (dot(L, N) <= 0.0)
                return false;

            rScattered = createRay(rec.pos + N * epsilon, normalize(L));

            float NoH = max(dot(N, H), 0.001);
            float NoV = max(dot(N, V), 0.001);
            float NoL = max(dot(N, L), 0.001);

            vec3 F = fresnelSchlick(max(dot(H, V), 0.0), rec.material.specColor);
            float D = D_GGX(NoH, roughness);
            float G = G_Smith(NoV, NoL, roughness);

            vec3 numerator = F * D * G;
            float denominator = 4.0 * max(dot(N, V), 0.001) * max(dot(N, L), 0.001);
            atten = numerator / max(denominator, 0.001);
        }
        else // Difuso
        {
            vec3 scatterDirection = N + randomUnitVector(gSeed);
            if (length(scatterDirection) < 1e-8)
                scatterDirection = N;

            rScattered = createRay(rec.pos + N * epsilon, normalize(scatterDirection));
            atten = rec.material.albedo * max(dot(N, rScattered.d), 0.0) / 3.141592;
        }

        return true;
    }


    if (rec.material.type == MT_DIELECTRIC)
    {
        vec3 outwardNormal;
        float niOverNt;
        float cosine;
        vec3 rdNorm = normalize(rIn.d);
        float dist = rec.t;
        bool isInside = dot(rdNorm, rec.normal) > 0.0;

        if (isInside)
        {
            outwardNormal = -rec.normal;
            niOverNt = rec.material.refIdx;
            cosine = dot(rdNorm, rec.normal);
        }
        else
        {
            outwardNormal = rec.normal;
            niOverNt = 1.0 / rec.material.refIdx;
            cosine = -dot(rdNorm, rec.normal);
        }

        if (isInside)
            atten = exp(-rec.material.refractColor * dist);
        else
            atten = vec3(1.0);

        vec3 refracted = refract(rdNorm, outwardNormal, niOverNt);
        bool canRefract = length(refracted) > 0.0001;

        float cosThetaT = dot(-normalize(refracted), outwardNormal);
        float reflectProb = canRefract ? schlick(cosThetaT, rec.material.refIdx) : 1.0;

        if (hash1(gSeed) < reflectProb)
        {
            vec3 reflected = reflect(rdNorm, rec.normal);
            vec3 scatteredDir = normalize(reflected + rec.material.roughness * randomInUnitSphere(gSeed));
            rScattered = createRay(rec.pos + rec.normal * epsilon, scatteredDir);
        }
        else
        {
            vec3 fuzzedRefracted = normalize(refracted + rec.material.roughness * randomInUnitSphere(gSeed));
            rScattered = createRay(rec.pos - outwardNormal * epsilon, fuzzedRefracted);
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



// Função que verifica a interseção entre um raio e uma esfera.
// Se houver interseção válida no intervalo [tMin, tMax], atualiza o HitRecord e retorna true.
// Caso contrário, retorna false.
bool hit_sphere(Sphere sph, Ray ray, float tMin, float tMax, inout HitRecord hit)
{
    // Vetor do centro da esfera até à origem do raio
    vec3 relOrigin = ray.o - sph.center;

    // Termo linear da equação quadrática (sem coeficiente a, pois o raio está normalizado)
    float a = dot(relOrigin, ray.d);

    // Discriminante da equação quadrática (sem coeficiente a²)
    float discriminant = a * a - (dot(relOrigin, relOrigin) - sph.radius * sph.radius);

    // Se o discriminante for negativo, não há interseção
    if (discriminant < 0.0) return false;

    // Verifica se o raio está fora da esfera inicialmente
    bool outside = dot(relOrigin, relOrigin) - sph.radius * sph.radius > 0.0;

    // Raiz quadrada do discriminante para cálculo da solução da equação
    float sqrtDisc = sqrt(discriminant);

    // Escolhe a solução correta da equação (mais próxima, dependendo se está fora ou dentro)
    float tCandidate = -a + (outside ? -sqrtDisc : sqrtDisc);

    // Verifica se o ponto de interseção está dentro dos limites válidos
    if (tCandidate < tMin || tCandidate > tMax)
        return false;

    // Calcula a posição e a normal no ponto de interseção
    vec3 intersection = ray.o + tCandidate * ray.d;
    vec3 normalVec = normalize(intersection - sph.center);

    // Inverte a normal se o raio intersecta o interior de uma esfera negativa
    if (sph.radius < 0.0)
        normalVec = -normalVec;

    // Atualiza os dados do HitRecord
    hit.t = tCandidate;
    hit.pos = intersection;
    hit.normal = normalVec;

    return true;
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