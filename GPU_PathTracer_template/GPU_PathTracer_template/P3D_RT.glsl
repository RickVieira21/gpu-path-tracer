/**
* ver hash functions em
* https://www.shadertoy.com/view/XlGcRh hash functions GPU
* http://www.jcgt.org/published/0009/03/02/
 */

 #include "./common.glsl"
 #iChannel0 "self"
 
 #define SCENE 0

bool hit_world(Ray r, float tmin, float tmax, inout HitRecord rec)
{
    bool hit = false;
    rec.t = tmax;

    #if SCENE == 0       //Shirley Weekend scene

        if(hit_quad(createQuad(vec3(-10.0, -0.05, 10.0), vec3(10.0, -0.05, 10.0), vec3(10.0, -0.05, -10.0), vec3(-10.0, -0.05, -10.0)), r, tmin, rec.t, rec))
        {
            hit = true;
            rec.material = createDiffuseMaterial(vec3(0.2));
        }

        if(hit_sphere(createSphere(vec3(-4.0, 1.0, 0.0), 1.0), r, tmin, rec.t, rec))
        {
            hit = true;
            rec.material = createDiffuseMaterial(vec3(0.2, 0.95, 0.1));
            //rec.material = createDiffuseMaterial(vec3(0.4, 0.2, 0.1));
        }

        if(hit_sphere(createSphere(vec3(4.0, 1.0, 0.0), 1.0),r,tmin,rec.t,rec))
        {
            hit = true;
            rec.material = createMetalMaterial(vec3(0.7, 0.6, 0.5), 0.0);
        }

        if(hit_sphere(createSphere(vec3(-1.5, 1.0, 0.0), 1.0),r,tmin,rec.t,rec))
        {
            hit = true;
            rec.material = createDielectricMaterial(vec3(0.0), 1.33, 0.0);
        }

        if(hit_sphere(createSphere(vec3(-1.5, 1.0, 0.0), -0.5),r,tmin,rec.t,rec))
        {
            hit = true;
            rec.material = createDielectricMaterial(vec3(0.0), 1.33, 0.0);
        }

        if(hit_sphere(createSphere(vec3(1.5, 1.0, 0.0), 1.0),r,tmin,rec.t,rec))
        {
            hit = true;
            rec.material = createDielectricMaterial(vec3(0.0, 0.9, 0.9), 1.5, 0.0);
        }
            
        int numxy = 5;
        
        for(int x = -numxy; x < numxy; ++x)
        {
            for(int y = -numxy; y < numxy; ++y)
            {
                float fx = float(x);
                float fy = float(y);
                float seed = fx + fy / 1000.0;
                vec3 rand1 = hash3(seed);
                vec3 center = vec3(fx + 0.9 * rand1.x, 0.2, fy + 0.9 * rand1.y);
                float chooseMaterial = rand1.z;
                if(distance(center, vec3(4.0, 0.2, 0.0)) > 0.9)
                {
                    if(chooseMaterial < 0.3)
                    {
                        vec3 center1 = center + vec3(0.0, hash1(gSeed) * 0.5, 0.0);
                        // diffuse
                        if(hit_movingSphere(createMovingSphere(center, center1, 0.2, 0.0, 1.0),r,tmin,rec.t,rec))
                        {
                            hit = true;
                            rec.material = createDiffuseMaterial(hash3(seed) * hash3(seed));
                        }
                    }
                    else if(chooseMaterial < 0.5)
                    {
                        // diffuse
                        if(hit_sphere(createSphere(center, 0.2),r,tmin,rec.t,rec))
                        {
                            hit = true;
                            rec.material = createDiffuseMaterial(hash3(seed) * hash3(seed));
                        }
                    }
                    else if(chooseMaterial < 0.7)
                    {
                        // metal
                        if(hit_sphere(createSphere(center, 0.2),r,tmin,rec.t,rec))
                        {
                            hit = true;
                            rec.material = createMetalMaterial((hash3(seed) + 1.0) * 0.5, 0.0);
                        }
                    }
                    else if(chooseMaterial < 0.9)
                    {
                        // metal
                        if(hit_sphere(createSphere(center, 0.2),r,tmin,rec.t,rec))
                        {
                            hit = true;
                            rec.material = createMetalMaterial((hash3(seed) + 1.0) * 0.5, hash1(seed));
                        }
                    }
                    else
                    {
                        // glass (Dielectric)
                        if(hit_sphere(createSphere(center, 0.2),r,tmin,rec.t,rec))
                        {
                            hit = true;
                            rec.material = createDielectricMaterial(hash3(seed), 1.33, 0.0);
                        }
                    }
                }
            }
        }
    #elif SCENE == 1 //from https://blog.demofox.org/2020/06/14/casual-shadertoy-path-tracing-3-fresnel-rough-refraction-absorption-orbit-camera/

        // diffuse floor
        
            vec3 A = vec3(-25.0f, -12.5f, 10.0f);
            vec3 B = vec3( 25.0f, -12.5f, 10.0f);
            vec3 C = vec3( 25.0f, -12.5f, -5.0f);
            vec3 D = vec3(-25.0f, -12.5f, -5.0f);

            if(hit_quad(createQuad(A, B, C, D), r, tmin, rec.t, rec))
            {
                hit = true;
                rec.material = createDiffuseMaterial(vec3(0.7));
            }

        //stripped background
        {
            vec3 A = vec3(-25.0f, -10.5f, -5.0f);
            vec3 B = vec3( 25.0f, -10.5f, -5.0f);
            vec3 C = vec3( 25.0f, -1.5f, -5.0f);
            vec3 D = vec3(-25.0f, -1.5f, -5.0f);
        
            if(hit_quad(createQuad(A, B, C, D), r, tmin, rec.t, rec))
            {
                hit = true;
                float shade = floor(mod(rec.pos.x, 1.0f) * 2.0f);
                rec.material = createDiffuseMaterial(vec3(shade));
            }
        }

        // ceiling piece above light
        
        {
            vec3 A = vec3(-7.5f, 12.5f, 5.0f);
            vec3 B = vec3( 7.5f, 12.5f, 5.0f);
            vec3 C = vec3( 7.5f, 12.5f, -5.0f);
            vec3 D = vec3(-7.5f, 12.5f, -5.0f);

            if(hit_quad(createQuad(A, B, C, D), r, tmin, rec.t, rec))
            {
                hit = true;
                rec.material = createDiffuseMaterial(vec3(0.7));
            }
        }    
       
        // light
        
        {
            vec3 A = vec3(-5.0f, 12.3f,  2.5f);
            vec3 B = vec3( 5.0f, 12.3f,  2.5f);
            vec3 C = vec3( 5.0f, 12.3f,  -2.5f);
            vec3 D = vec3(-5.0f, 12.3f,  -2.5f);

             if(hit_quad(createQuad(A, B, C, D), r, tmin, rec.t, rec))
            {
                hit = true;
                rec.material = createDiffuseMaterial(vec3(0.0));
                rec.material.emissive = vec3(1.0f, 0.9f, 0.9f) * 20.0f;
            }
        }
 
        const int c_numSpheres = 7;
        for (int sphereIndex = 0; sphereIndex < c_numSpheres; ++sphereIndex)
        {
            vec3 center = vec3(-18.0 + 6.0 * float(sphereIndex), -8.0, 0.0);
            if(hit_sphere(createSphere(center, 2.8),r,tmin,rec.t,rec))
            {
                hit = true;
                float r = float(sphereIndex) / float(c_numSpheres-1) * 0.1f;
                rec.material = createDielectricMaterial(vec3(0.0, 0.5, 1.0), 1.1, r);
            }
        }

    #elif SCENE == 2
    #elif SCENE == 3
    #endif

    return hit;
}

vec3 directlighting(pointLight pl, Ray r, HitRecord rec) {
    float shininess = 1.0;
    vec3 normal = normalize(rec.normal);
    vec3 surfacePos = rec.pos + normal * epsilon;

    vec3 toLight = pl.pos - surfacePos;
    float lightDist = length(toLight);
    vec3 lightDir = normalize(toLight);

    // Shadow ray
    Ray shadowRay = createRay(surfacePos, lightDir);
    HitRecord shadowRec;

    // Check if something blocks the light
    if (hit_world(shadowRay, 0.001, lightDist - 0.001, shadowRec)) {
        return vec3(0.0); // in shadow → no direct lighting
    }

    // If visible, compute light
    vec3 diffuse = rec.material.albedo * max(dot(normal, lightDir), 0.0);
    vec3 halfway = normalize(lightDir - r.d);
    vec3 specular = pl.color * rec.material.emissive * pow(max(dot(halfway, normal), 0.0), shininess);

    return diffuse + specular;
}


#define SHADOW_SAMPLES 4

vec3 samplePointOnQuad(vec3 A, vec3 B, vec3 C, vec3 D, vec2 uv) {
    // Bilinear interpolation on quad surface
    vec3 AB = mix(A, B, uv.x);
    vec3 DC = mix(D, C, uv.x);
    return mix(AB, DC, uv.y);
}

vec3 directLightingEmissiveQuad(vec3 hitPoint, vec3 viewDir, vec3 normal, Material mat)
{
    // Only skip direct lighting *on* emissive materials themselves, but allow indirect lighting on others
    if(length(mat.emissive) > 0.001)
        return vec3(0.0);

    vec3 A = vec3(-5.0, 12.3, 2.5);
    vec3 B = vec3(5.0, 12.3, 2.5);
    vec3 C = vec3(5.0, 12.3, -2.5);
    vec3 D = vec3(-5.0, 12.3, -2.5);
    vec3 lightEmission = vec3(1.0, 0.9, 0.9) * 2.0;

    vec3 colorAcc = vec3(0.0);
    int gridSize = 2; // 2x2 samples for soft shadows

    for (int s = 0; s < gridSize * gridSize; ++s) {
        int i = s / gridSize;
        int j = s % gridSize;

        vec2 uv = (vec2(float(i), float(j)) + 0.5) / float(gridSize);
        vec3 samplePos = samplePointOnQuad(A, B, C, D, uv);

        vec3 lightVec = samplePos - hitPoint;
        float distToLight = length(lightVec);
        vec3 lightDir = normalize(lightVec);

        vec3 offsetHitPoint = hitPoint + normal * 0.001;

        Ray shadowRay = createRay(offsetHitPoint, lightDir);
        HitRecord shadowRec;

        float maxShadowDist = max(distToLight - 0.001, 0.001);
        bool inShadow = hit_world(shadowRay, 0.001, maxShadowDist, shadowRec);

        if (!inShadow) {
            float NdotL = max(dot(normal, lightDir), 0.0);

            vec3 halfway = normalize(lightDir + viewDir);
            float NdotH = max(dot(normal, halfway), 0.0);

            // Map roughness to shininess exponent (simple mapping)
            float shininess = 1.0 / max(mat.roughness, 0.001);
            float specularFactor = pow(NdotH, shininess);

            vec3 diffuse = mat.albedo * lightEmission * NdotL;
            vec3 specular = mat.specColor * lightEmission * specularFactor;

            colorAcc += diffuse + specular;
        }
    }

    return colorAcc / float(gridSize * gridSize);
}




#define MAX_BOUNCES 10

vec3 rayColor(Ray r)
{
    HitRecord rec;
    vec3 col = vec3(0.0);
    vec3 throughput = vec3(1.0);

    for (int i = 0; i < MAX_BOUNCES; ++i)
    {
        if (hit_world(r, 0.001, 10000.0, rec))
        {
            vec3 viewDir = normalize(-r.d);

            //col += directlighting(createPointLight(vec3(-10.0, 15.0, 0.0), vec3(1.0)), r, rec) * throughput;
            //col += directlighting(createPointLight(vec3(8.0, 15.0, 3.0), vec3(1.0)), r, rec) * throughput;
            //col += directlighting(createPointLight(vec3(1.0, 15.0, -9.0), vec3(1.0)), r, rec) * throughput;

            vec3 lighting = vec3(0.0);
            // Add emissive quad lighting with soft shadows
            lighting += directLightingEmissiveQuad(rec.pos, viewDir, rec.normal, rec.material);

            col += throughput * lighting;

            // Add emission from the material itself
            col += throughput * rec.material.emissive;

            Ray scatterRay;
            vec3 attenuation;
            if (scatter(r, rec, attenuation, scatterRay))
            {
                throughput *= attenuation;
                if(any(lessThan(throughput, vec3(0.0)))) break; // prevent negative throughput
                r = scatterRay;
            }
            else
            {
                break;
            }
        }
        else
        {
            float t = 0.8 * (r.d.y + 1.0);
            col += throughput * mix(vec3(1.0), vec3(0.5, 0.7, 1.0), t);
            break;
        }
    }

    return col;
}


#define PI 3.14159265358979
#define MAX_SAMPLES 10000.0

void main()
{
    gSeed = float(baseHash(floatBitsToUint(gl_FragCoord.xy))) / float(0xffffffffU) + iTime;

   //Orbital Cam
    vec2 mouse = iMouse.xy / iResolution.xy;
    mouse = clamp(mouse, 0.001, 0.999); // evitar extremos 0 ou 1

    float radius = 10.0 + 6.0 * (1.0 - mouse.y);           // zoom controlado no eixo Y
    float alpha = mouse.x * 2.0 * PI;                     // ângulo horizontal (órbita)
    float beta = mix(-PI * 0.25, PI * 0.25, mouse.y);     // ângulo vertical limitado

    vec3 camTarget = vec3(0.0); // centro da cena

    vec3 camPos;
    camPos.x = camTarget.x + radius * cos(beta) * sin(alpha);
    camPos.y = camTarget.y + radius * sin(beta);
    camPos.z = camTarget.z + radius * cos(beta) * cos(alpha);


    float fovy = 60.0;
    float aperture = 0.0;
    float distToFocus = 1.0;
    float time0 = 0.0;
    float time1 = 1.0;
    Camera cam = createCamera(
        camPos,
        camTarget,
        vec3(0.0, 1.0, 0.0),    // world up vector
        fovy,
        iResolution.x / iResolution.y,
        aperture,
        distToFocus,
        time0,
        time1);

    //usa-se o 4 canal de cor para guardar o numero de samples e não o iFrame pois quando se mexe o rato faz-se reset

    vec4 prev = texture(iChannel0, gl_FragCoord.xy / iResolution.xy);
    vec3 prevLinear = toLinear(prev.xyz);  

    vec2 ps = gl_FragCoord.xy + hash2(gSeed);
    //vec2 ps = gl_FragCoord.xy;
    vec3 color = rayColor(getRay(cam, ps));

    if(iMouseButton.x != 0.0 || iMouseButton.y != 0.0)
    {
        gl_FragColor = vec4(toGamma(color), 1.0);  //samples number reset = 1
        return;
    }
    if(prev.w > MAX_SAMPLES)   
    {
        gl_FragColor = prev;
        return;
    }

    float w = prev.w + 1.0;
    color = mix(prevLinear, color, 1.0/w);
    gl_FragColor = vec4(toGamma(color), w);
}