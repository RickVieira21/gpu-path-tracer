/**
 * common.glsl
 * Common types and functions used for ray tracing.
 */

const float pi = 3.14159265358979;
const float epsilon = 0.001;
const float RECIPROCAL_PI = 1.0 / 3.141592653589793;

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

// Fresnel-Schlick: aproximação eficiente da equação de Fresnel
// cosTheta: cosseno do ângulo entre a normal e o vetor de visão (V.H)
// F0: refletância no ângulo normal de incidência (0º), tipicamente 0.04 para dielétricos, ou a cor metálica para metais
vec3 fresnelSchlick(float cosTheta, vec3 F0) {
    return F0 + (1.0 - F0) * pow(1.0 - cosTheta, 5.0);
}


// D_GGX: função de distribuição normal (NDF) para o modelo GGX 
// NoH: cosseno entre a normal da superfície (N) e o vetor H (half-way entre L e V)
// roughness: rugosidade do material (entre 0 e 1)
float D_GGX(float NoH, float roughness) {
    float alpha = roughness * roughness;    // converte rugosidade perceptual para "alpha"
    float alpha2 = alpha * alpha;           // alpha ao quadrado
    float NoH2 = NoH * NoH;                 // cos^2 entre N e H
    float b = (NoH2 * (alpha2 - 1.0) + 1.0); // denom. comum
    return alpha2 / (pi * b * b);           // distribuição GGX normalizada
}


// G1_GGX_Schlick: aproximação da função de geometria de visibilidade para um dos lados (L ou V)
// NoV: cosseno entre N e o vetor de visão ou luz
float G1_GGX_Schlick(float NoV, float roughness) {
    float r = 0.5 + 0.5 * roughness;
    float k = (r * r) / 2.0;                // fator de suavização
    float denom = NoV * (1.0 - k) + k;      // denom. de Schlick
    return max(NoV, 0.001) / denom;         // evita divisões por zero
}


// G_Smith: função de geometria completa (Smith) para GGX
// Combina os dois lados (visão e luz) com a aproximação de Schlick
float G_Smith(float NoV, float NoL, float roughness) {
    float g1_l = G1_GGX_Schlick(NoL, roughness); // lado da luz
    float g1_v = G1_GGX_Schlick(NoV, roughness); // lado da visão
    return g1_l * g1_v; // visibilidade combinada
}


// sampleGGX: gera um vetor H amostrado segundo a distribuição GGX
// (ou seja vetor médio entre o vetor de visão V (direção da câmara) e o vetor da luz L)
// V: vetor de visão (normalizado, apontando para a câmara)
vec3 sampleGGX(vec3 V, float roughness, inout float seed)
{
    float a = roughness * roughness; // alpha (rugosidade convertida)

    // Gera coordenadas esféricas (phi, theta) segundo GGX
    float phi = 2.0 * pi * hash1(seed); // ângulo aleatório
    float cosTheta = sqrt((1.0 - hash1(seed)) / (1.0 + (a*a - 1.0) * hash1(seed))); // GGX theta
    float sinTheta = sqrt(1.0 - cosTheta * cosTheta); // sen de theta 

    // Converte coordenadas esféricas para vetor no espaço local (tangente ao V)
    vec3 H = vec3(cos(phi) * sinTheta, sin(phi) * sinTheta, cosTheta);

    // Cria uma base ortonormal em torno de V para transformar H
    // Isto permite amostrar H no hemisfério orientado pela direção V
    vec3 up = abs(V.z) < 0.999 ? vec3(0.0, 0.0, 1.0) : vec3(1.0, 0.0, 0.0); // evita degenerescência
    vec3 tangentX = normalize(cross(up, V));   // primeiro vetor tangente
    vec3 tangentY = cross(V, tangentX);        // segundo vetor tangente ortogonal

    // Converte H do espaço tangente para o espaço do mundo (orientado por V)
    return normalize(H.x * tangentX + H.y * tangentY + H.z * V);
}


vec3 brdfMicrofacet(in vec3 L, in vec3 V, in vec3 N,
                    in float metallic, in float roughness,
                    in vec3 baseColor, in float reflectance)
{
    // Halfway vector entre V e L
    vec3 H = normalize(V + L);

    // Produtos escalares relevantes para a BRDF
    float NoV = clamp(dot(N, V), 0.0, 1.0);  // cos(θ_view)
    float NoL = clamp(dot(N, L), 0.0, 1.0);  // cos(θ_light)
    float NoH = clamp(dot(N, H), 0.0, 1.0);  // cos(θ_half)
    float VoH = clamp(dot(V, H), 0.0, 1.0);  // cos(θ_fresnel)

    // F0: refletância na incidência normal (usada no fresnel)
    vec3 f0 = vec3(0.16) * (reflectance * reflectance);  // ex: 0.04 para plásticos
    f0 = mix(f0, baseColor, metallic); // metais usam baseColor como F0

    // Termos da BRDF especular de Cook-Torrance
    vec3 F = fresnelSchlick(VoH, f0);         // Fresnel (reflexão angular)
    float D = D_GGX(NoH, roughness);          // Distribuição de microfacetos (GGX)
    float G = G_Smith(NoV, NoL, roughness);   // Termo de geometria (oclusão)

    // Componente especular final
    vec3 spec = (F * D * G) / (4.0 * max(NoV, 0.001) * max(NoL, 0.001));

    // Componente difusa (zero para metais)
    vec3 rhoD = baseColor * (1.0 - metallic); // sem difuso para metal
    vec3 diff = rhoD * RECIPROCAL_PI;         // lambertiano: 1/π

    // BRDF = difusa + especular
    return diff + spec;
}



// Old scatter: calcula a dispersão do raio com base no tipo de material
bool scatter(Ray rIn, HitRecord rec, out vec3 atten, out Ray rScattered)
{
    // --- Difuso Lambertiano ---
    if (rec.material.type == MT_DIFFUSE)
    {
        // Calcula uma direção de dispersão aleatória a partir da normal (Lambertian)
        vec3 scatterDirection = rec.normal + randomUnitVector(gSeed);

        // Evita direções próximas de zero (artefatos)
        if (length(scatterDirection) < 1e-8)
            scatterDirection = rec.normal;

        // Cria o novo raio disperso a partir da posição da interseção
        rScattered = createRay(rec.pos, scatterDirection);

        // Atenuação da energia com base no albedo e no ângulo entre o raio e a normal
        atten = rec.material.albedo * max(dot(normalize(rScattered.d), rec.normal), 0.0) / 3.141592;

        return true; // O raio foi disperso
    }

    // --- Metal (reflexão com rugosidade) ---
    if (rec.material.type == MT_METAL)
    {
        // Calcula o vetor refletido
        vec3 reflected = reflect(normalize(rIn.d), rec.normal);

        // Adiciona ruído baseado na rugosidade → cria "fuzziness"
        vec3 scatteredDir = reflected + rec.material.roughness * randomInUnitSphere(gSeed);

        // Cria o raio refletido a partir da posição da interseção
        rScattered = createRay(rec.pos, scatteredDir);

        // Atenuação com a cor especular do material
        atten = rec.material.specColor;

        // Só reflete se estiver a sair da superfície (dot positivo com a normal)
        return dot(rScattered.d, rec.normal) > 0.0;
    }

    // --- Dielétrico ---
    if (rec.material.type == MT_DIELECTRIC)
    {
        vec3 outwardNormal;
        float niOverNt;
        float cosine;
        vec3 rdNorm = normalize(rIn.d); // Direção normalizada do raio de entrada
        float dist = rec.t; // Distância percorrida até ao ponto de interseção

        // Verifica se o raio está dentro ou fora do meio
        bool isInside = dot(rdNorm, rec.normal) > 0.0;

        if (isInside)
        {
            // O raio está a sair do meio → invertendo a normal
            outwardNormal = -rec.normal;
            niOverNt = rec.material.refIdx; // Índice de refração do meio

            cosine = dot(rdNorm, rec.normal); // Coseno do ângulo de incidência
        }
        else
        {
            // O raio está a entrar no meio → normal permanece
            outwardNormal = rec.normal;
            niOverNt = 1.0 / rec.material.refIdx; // Ar para meio refrativo

            cosine = -dot(rdNorm, rec.normal);
        }

        // Aplica a lei de Beer para absorção se estiver dentro do meio
        if (isInside)
            atten = exp(-rec.material.refractColor * dist); // Atenuação exponencial
        else
            atten = vec3(1.0); // Sem atenuação ao entrar no meio

        // Tenta calcular o vetor refratado
        vec3 refracted = refract(rdNorm, outwardNormal, niOverNt);
        bool canRefract = length(refracted) > 0.0001; // Confirma que refratou

        // Coseno do ângulo da refratada com a normal externa (para Schlick)
        float cosThetaT = dot(-normalize(refracted), outwardNormal); 

        // Usa aproximação de Schlick para refletividade
        float reflectProb = canRefract ? schlick(cosThetaT, rec.material.refIdx) : 1.0;

        // Escolhe entre refletir ou refratar com base numa amostra aleatória
        if (hash1(gSeed) < reflectProb)
        {
            // Reflete com rugosidade
            vec3 reflected = reflect(rdNorm, rec.normal);
            vec3 scatteredDir = normalize(reflected + rec.material.roughness * randomInUnitSphere(gSeed));

            // Cria o raio refletido ligeiramente afastado da superfície
            rScattered = createRay(rec.pos + rec.normal * epsilon, scatteredDir);
        }
        else
        {
            // Refrata com alguma rugosidade (fuzzing)
            vec3 fuzzedRefracted = normalize(refracted + rec.material.roughness * randomInUnitSphere(gSeed));

            // Cria o raio refratado ligeiramente para dentro
            rScattered = createRay(rec.pos - outwardNormal * epsilon, fuzzedRefracted);
        }

        return true; // Raio foi disperso (refletido ou refratado)
    }

    return false;
}


//New Scatter 
/*
bool scatter(Ray rIn, HitRecord rec, out vec3 atten, out Ray rScattered)
{
    // ======== DIFUSO (Lambertiano) ========
    if (rec.material.type == MT_DIFFUSE)
    {
        // Direção de dispersão baseada no modelo Lambertiano: normal + vetor aleatório unitário
        vec3 scatterDirection = rec.normal + randomUnitVector(gSeed);

        // Evita vetores degenerados (muito próximos de zero)
        if (length(scatterDirection) < 1e-8)
            scatterDirection = rec.normal;

        // Raio dispersado é originado na superfície e segue na direção escolhida
        rScattered = createRay(rec.pos, scatterDirection);

        // Atenuação baseada na reflectância difusa (albedo), normalizada por PI (modelo Lambertiano)
        atten = rec.material.albedo * max(dot(normalize(rScattered.d), rec.normal), 0.0) / pi;
        return true;
    }

    // ======== METAL (Reflexão com microfacet) ========
    if (rec.material.type == MT_METAL)
    {
        vec3 V = normalize(-rIn.d);         // Vetor de visualização (do ponto para a câmara)
        vec3 N = rec.normal;
        float roughness = rec.material.roughness;
        vec3 specColor = rec.material.specColor;

        // H é o vetor halfway (entre V e L), amostrado com distribuição GGX para rugosidade
        vec3 H = sampleGGX(V, roughness, gSeed);

        // Calcula L refletindo V em relação a H
        vec3 L = reflect(-V, H); 

        // Se L estiver do lado errado da superfície, a reflexão não ocorre
        if (dot(L, N) <= 0.0)
            return false;

        // Cria o raio refletido com pequeno offset para evitar self-intersection
        rScattered = createRay(rec.pos + N * epsilon, normalize(L));

        // Usa BRDF de microfacet com metal (metallic = 1.0, baseColor = specColor)
        atten = brdfMicrofacet(L, V, N, 1.0, roughness, specColor, 1.0) * max(dot(N, L), 0.0);

        return true;
    }

    // ======== PLÁSTICO (Reflexão + Difuso com Fresnel) ========
    if (rec.material.type == MT_PLASTIC)
    {
        vec3 V = normalize(-rIn.d);
        vec3 N = rec.normal;
        float roughness = rec.material.roughness;
        vec3 baseColor = rec.material.albedo;
        float reflectance = rec.material.refIdx; 

        // Cálculo de Fresnel com base na cor especular
        float cosTheta = max(dot(N, V), 0.0);
        float fresnel = fresnelSchlick(cosTheta, rec.material.specColor).r;

        if (hash1(gSeed) < fresnel)
        {
            // Reflexão especular com microfacet
            vec3 H = sampleGGX(N, roughness, gSeed);
            vec3 L = reflect(-V, H);

            if (dot(L, N) <= 0.0)
                return false;

            rScattered = createRay(rec.pos + N * epsilon, normalize(L));

            // Reflexão especular com BRDF microfacet (metallic = 0.0)
            atten = brdfMicrofacet(L, V, N, 0.0, roughness, baseColor, reflectance) * max(dot(N, L), 0.0);
        }
        else
        {
            // Reflexão difusa (como Lambert) com randomização
            vec3 scatterDirection = N + randomUnitVector(gSeed);
            if (length(scatterDirection) < 1e-8)
                scatterDirection = N;

            rScattered = createRay(rec.pos + N * epsilon, normalize(scatterDirection));
            vec3 L = rScattered.d;

            // Atenuação difusa tratada também com BRDF para consistência
            atten = brdfMicrofacet(L, V, N, 0.0, roughness, baseColor, reflectance) * max(dot(N, L), 0.0);
        }

        return true;
    }

    // ======== DIELÉTRICO (Transparente: Refração + Reflexão) ========
    if (rec.material.type == MT_DIELECTRIC)
    {
        vec3 outwardNormal;
        float niOverNt;     // Razão de índices de refração
        float cosine;
        vec3 rdNorm = normalize(rIn.d);
        float dist = rec.t;

        // Verifica se o raio está dentro do material
        bool isInside = dot(rdNorm, rec.normal) > 0.0;

        if (isInside)
        {
            outwardNormal = -rec.normal;
            niOverNt = rec.material.refIdx;           // de dentro para fora
            cosine = dot(rdNorm, rec.normal);
        }
        else
        {
            outwardNormal = rec.normal;
            niOverNt = 1.0 / rec.material.refIdx;      // de fora para dentro
            cosine = -dot(rdNorm, rec.normal);
        }

        // Atenuação por absorção (apenas quando o raio atravessa o meio)
        if (isInside)
            atten = exp(-rec.material.refractColor * dist); // absorção volumétrica
        else
            atten = vec3(1.0); // sem absorção fora do meio

        vec3 refracted = refract(rdNorm, outwardNormal, niOverNt); // tenta calcular o vetor refratado com base na normal e razão dos índices de refração
        bool canRefract = length(refracted) > 0.0001;              // verifica se foi possível refratar (comprimento próximo de zero indica falha)

        float cosThetaT = dot(-normalize(refracted), outwardNormal); // calcula o cosseno do ângulo entre o vetor refratado (invertido) e a normal
        float reflectProb = canRefract ? schlick(cosThetaT, rec.material.refIdx) : 1.0; // usa a aproximação de Schlick para estimar a proporção refletida


        if (hash1(gSeed) < reflectProb)
        {
            // Reflexão com fuzziness (se houver rugosidade)
            vec3 reflected = reflect(rdNorm, rec.normal);
            vec3 scatteredDir = normalize(reflected + rec.material.roughness * randomInUnitSphere(gSeed));
            rScattered = createRay(rec.pos + rec.normal * epsilon, scatteredDir);
        }
        else
        {
            // Refração com fuzziness
            vec3 fuzzedRefracted = normalize(refracted + rec.material.roughness * randomInUnitSphere(gSeed));
            rScattered = createRay(rec.pos - outwardNormal * epsilon, fuzzedRefracted);
        }

        return true;
    }

    // Caso o tipo de material não seja reconhecido
    return false;
}
*/





struct Triangle {vec3 a; vec3 b; vec3 c; };

Triangle createTriangle(vec3 v0, vec3 v1, vec3 v2)
{
    Triangle t;
    t.a = v0; t.b = v1; t.c = v2;
    return t;
}


bool hit_triangle(Triangle t, Ray r, float tmin, float tmax, out HitRecord rec)
{
    vec3 v0 = t.a;
    vec3 v1 = t.b;
    vec3 v2 = t.c;

    // Calcula as arestas do triângulo
    vec3 edge1 = v1 - v0;
    vec3 edge2 = v2 - v0;

    // Produto vetorial entre a direção do raio e a segunda aresta
    vec3 h = cross(r.d, edge2);
    float a = dot(edge1, h); // Determinante do sistema

    // Comentado: verificação se o raio é paralelo ao triângulo
    // if (abs(a) < 1e-8) return false;

    float f = 1.0 / a;
    vec3 s = r.o - v0;
    float u = f * dot(s, h); // Coordenada baricêntrica u

    // Verifica se o ponto está fora do triângulo
    if (u < 0.0 || u > 1.0)
        return false;

    vec3 q = cross(s, edge1);
    float v = f * dot(r.d, q); // Coordenada baricêntrica v

    // Verifica se a interseção está fora do triângulo
    if (v < 0.0 || u + v > 1.0)
        return false;

    float t_intersect = f * dot(edge2, q); // Distância até à interseção

    // Verifica se está dentro do intervalo válido
    if (t_intersect < tmin || t_intersect > tmax)
        return false;

    // Interseção válida: preenche o HitRecord
    rec.t = t_intersect;
    rec.pos = r.o + t_intersect * r.d;
    rec.normal = normalize(cross(edge1, edge2)); // Normal do triângulo

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


vec3 center(MovingSphere mvsphere, float time)
{
    // Interpola linearmente a posição da esfera com base no tempo
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


bool hit_movingSphere(MovingSphere s, Ray r, float tmin, float tmax, out HitRecord rec)
{
    // 1. Interpola o centro da esfera com base no tempo do raio
    float timeFraction = (r.t - s.time0) / (s.time1 - s.time0);
    vec3 center = mix(s.center0, s.center1, clamp(timeFraction, 0.0, 1.0));

    // 2. Interseção raio-esfera (como no caso estático)
    vec3 oc = r.o - center;
    float a = dot(r.d, r.d);
    float b = dot(oc, r.d);
    float c = dot(oc, oc) - s.radius * s.radius;
    float discriminant = b * b - a * c;

    if (discriminant > 0.0) {
        float sqrtD = sqrt(discriminant);

        // Primeira raiz
        float temp = (-b - sqrtD) / a;
        if (temp < tmax && temp > tmin) {
            rec.t = temp;
            rec.pos = pointOnRay(r, rec.t);
            rec.normal = (rec.pos - center) / s.radius;

            // Ajusta normal se estiver a entrar na esfera
            bool frontFace = dot(r.d, rec.normal) < 0.0;
            rec.normal = frontFace ? rec.normal : -rec.normal;

            return true;
        }

        // Segunda raiz (caso a primeira não seja válida)
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

    return false; // Sem interseção válida
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