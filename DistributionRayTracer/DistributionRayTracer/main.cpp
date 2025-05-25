 ///////////////////////////////////////////////////////////////////////
//
// P3D Course
// (c) 2025 by João Madeiras Pereira
//Distribution Ray Tracing P3F scenes and drawing points with Modern OpenGL
// It explores parallelism through OMP
//
///////////////////////////////////////////////////////////////////////
#include <omp.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string.h>
#include <stdio.h>
#include <chrono>
#include <conio.h>

#include <GL/glew.h>
#include <GL/freeglut.h>
#include <IL/il.h>

#include "scene.h"
#include "rayAccelerator.h"
#include "maths.h"
#include "macros.h"

//Enable OpenGL drawing.  
bool drawModeEnabled = true;
bool P3F_scene = true; //choose between P3F scene or the built-in Peter Shirley scene
bool Progressive_flg = false;

#define MAX_DEPTH 4  //number of bounces

#define CAPTION "Accel Distribution RT"
#define VERTEX_COORD_ATTRIB 0
#define COLOR_ATTRIB 1
#define MAX_SAMPLES 10000

// Frame counting and FPS computation
unsigned int FramesPerSecond = 0, FPS = 0;
unsigned long FrameCount = 1;

// Current Camera Position
float camX, camY, camZ;

//Original Camera position;
Vector Eye;

// Mouse Tracking Variables
int startX, startY, tracking = 0;

// Camera Spherical Coordinates
float alpha = 0.0f;
float _beta = 0.0f;
float r = 4.0f;

// Color Gamma correction
double invGamma = 1.0f / 2.2f;

// Points defined by 2 attributes: positions which are stored in vertices array and colors which are stored in colors array
float *colors;
float *vertices;
int size_vertices;
int size_colors;

//Array of Pixels to be stored in a file by using DevIL library
uint8_t *img_Data;

GLfloat m[16];  //projection matrix initialized by ortho function

GLuint VaoId;
GLuint VboId[2];

GLuint VertexShaderId, FragmentShaderId, ProgramId;
GLint UniformId;

Scene* scene = NULL;
Grid* grid_ptr = NULL;
BVH* bvh_ptr = NULL;

int RES_X, RES_Y;


int WindowHandle = 0;

bool AA = false;
unsigned int spp; //samples per pixel
bool DOF = false;
bool SoftShadows = false;

accelerator Accel_Struct = NONE;

/////////////////////////////////////////////////////////////////////// ERRORS

bool isOpenGLError() {
	bool isError = false;
	GLenum errCode;
	const GLubyte *errString;
	while ((errCode = glGetError()) != GL_NO_ERROR) {
		isError = true;
		errString = gluErrorString(errCode);
		std::cerr << "OpenGL ERROR [" << errString << "]." << std::endl;
	}
	return isError;
}

void checkOpenGLError(std::string error)
{
	if(isOpenGLError()) {
		std::cerr << error << std::endl;
		exit(EXIT_FAILURE);
	}
}

/////////////////////////////////////////////////////////////////////// SHADERs

const GLchar* VertexShader =
{
	"#version 430 core\n"

	"in vec2 in_Position;\n"
	"in vec3 in_Color;\n"
	"uniform mat4 Matrix;\n"
	"out vec4 color;\n"

	"void main(void)\n"
	"{\n"
	"	vec4 position = vec4(in_Position, 0.0, 1.0);\n"
	"	color = vec4(in_Color, 1.0);\n"
	"	gl_Position = Matrix * position;\n"

	"}\n"
};

const GLchar* FragmentShader =
{
	"#version 430 core\n"

	"in vec4 color;\n"
	"out vec4 out_Color;\n"

	"void main(void)\n"
	"{\n"
	"	out_Color = color;\n"
	"}\n"
};

void createShaderProgram()
{
	VertexShaderId = glCreateShader(GL_VERTEX_SHADER);
	glShaderSource(VertexShaderId, 1, &VertexShader, 0);
	glCompileShader(VertexShaderId);

	FragmentShaderId = glCreateShader(GL_FRAGMENT_SHADER);
	glShaderSource(FragmentShaderId, 1, &FragmentShader, 0);
	glCompileShader(FragmentShaderId);

	ProgramId = glCreateProgram();
	glAttachShader(ProgramId, VertexShaderId);
	glAttachShader(ProgramId, FragmentShaderId);

	glBindAttribLocation(ProgramId, VERTEX_COORD_ATTRIB, "in_Position");
	glBindAttribLocation(ProgramId, COLOR_ATTRIB, "in_Color");
	
	glLinkProgram(ProgramId);
	UniformId = glGetUniformLocation(ProgramId, "Matrix");

	checkOpenGLError("ERROR: Could not create shaders.");
}

void destroyShaderProgram()
{
	glUseProgram(0);
	glDetachShader(ProgramId, VertexShaderId);
	glDetachShader(ProgramId, FragmentShaderId);

	glDeleteShader(FragmentShaderId);
	glDeleteShader(VertexShaderId);
	glDeleteProgram(ProgramId);

	checkOpenGLError("ERROR: Could not destroy shaders.");
}

/////////////////////////////////////////////////////////////////////// VAOs & VBOs


void createBufferObjects()
{
	glGenVertexArrays(1, &VaoId);
	glBindVertexArray(VaoId);
	glGenBuffers(2, VboId);
	glBindBuffer(GL_ARRAY_BUFFER, VboId[0]);

	/* Só se faz a alocação dos arrays glBufferData (NULL), e o envio dos pontos para a placa gráfica
	é feito na drawPoints com GlBufferSubData em tempo de execução pois os arrays são GL_DYNAMIC_DRAW */
	glBufferData(GL_ARRAY_BUFFER, size_vertices, NULL, GL_DYNAMIC_DRAW);
	glEnableVertexAttribArray(VERTEX_COORD_ATTRIB);
	glVertexAttribPointer(VERTEX_COORD_ATTRIB, 2, GL_FLOAT, 0, 0, 0);
	
	glBindBuffer(GL_ARRAY_BUFFER, VboId[1]);
	glBufferData(GL_ARRAY_BUFFER, size_colors, NULL, GL_DYNAMIC_DRAW);
	glEnableVertexAttribArray(COLOR_ATTRIB);
	glVertexAttribPointer(COLOR_ATTRIB, 3, GL_FLOAT, 0, 0, 0);
	
// unbind the VAO
	glBindVertexArray(0);
    glBindBuffer(GL_ARRAY_BUFFER, 0);
//	glDisableVertexAttribArray(VERTEX_COORD_ATTRIB); 
//	glDisableVertexAttribArray(COLOR_ATTRIB);
	checkOpenGLError("ERROR: Could not create VAOs and VBOs.");
}

void destroyBufferObjects()
{
	glDisableVertexAttribArray(VERTEX_COORD_ATTRIB);
	glDisableVertexAttribArray(COLOR_ATTRIB);
	
	glBindBuffer(GL_ARRAY_BUFFER, 0);
	glBindVertexArray(0);

	glDeleteBuffers(1, VboId);
	glDeleteVertexArrays(1, &VaoId);
	checkOpenGLError("ERROR: Could not destroy VAOs and VBOs.");
}

void drawPoints()
{
	glClear(GL_COLOR_BUFFER_BIT);
	
	glBindVertexArray(VaoId);
	glUseProgram(ProgramId);

	glBindBuffer(GL_ARRAY_BUFFER, VboId[0]);
	glBufferSubData(GL_ARRAY_BUFFER, 0, size_vertices, vertices);
	glBindBuffer(GL_ARRAY_BUFFER, VboId[1]);
	glBufferSubData(GL_ARRAY_BUFFER, 0, size_colors, colors);

	glUniformMatrix4fv(UniformId, 1, GL_FALSE, m);
	glDrawArrays(GL_POINTS, 0, RES_X*RES_Y);
	//glFinish();
	glUseProgram(0);
	glBindVertexArray(0);
	glBindBuffer(GL_ARRAY_BUFFER, 0);

	checkOpenGLError("ERROR: Could not draw scene.");
	
}

ILuint saveImgFile(const char *filename) {
	ILuint ImageId;

	ilEnable(IL_FILE_OVERWRITE);
	ilGenImages(1, &ImageId);
	ilBindImage(ImageId);

	ilTexImage(RES_X, RES_Y, 1, 3, IL_RGB, IL_UNSIGNED_BYTE, img_Data /*Texture*/);
	ilSaveImage(filename);

	ilDisable(IL_FILE_OVERWRITE);
	ilDeleteImages(1, &ImageId);
	if (ilGetError() != IL_NO_ERROR)return ilGetError();

	return IL_NO_ERROR;
}

/////////////////////////////////////////////////////////////////////// CALLBACKS

void timer(int value)
{
	FramesPerSecond = FPS;
	if(!Progressive_flg) {
		std::ostringstream oss;
		oss << CAPTION << ": " << FramesPerSecond << " FPS @ (" << RES_X << "x" << RES_Y << ")";
		std::string s = oss.str();
		glutSetWindow(WindowHandle);
		glutSetWindowTitle(s.c_str());
	}
	FPS = 0;
	glutTimerFunc(1000, timer, 0);
}


///////////////////////////////////////////////////YOUR CODE HERE////////////////////////////////////////////////////////////////////////

Vector random_in_unit_sphere() {
	Vector p;
	do {
		p = Vector(rand_float() * 2.0f - 1.0f,
			rand_float() * 2.0f - 1.0f,
			rand_float() * 2.0f - 1.0f);
	} while (p.lengthSquared() >= 1.0f);
	return p;
}


Color rayTracing(Ray ray, int depth, float ior_1, Vector lightSample)  //index of refraction of medium 1 where the ray is travelling
{
	Color color_Acc; //Class constructor init the color with zero

	Object* hitObj = NULL; //nearest object
	HitRecord closestHit;  //isHit=false and t=FLT_MAX
	Vector hitPoint; //closest hit point
	Vector N;
	HitRecord auxRec;
	bool skybox_flg;
	int num_lights = scene->getNumLights();
	skybox_flg = scene->GetSkyBoxFlg();
	Accel_Struct = scene->GetAccelStruct();   //Type of acceleration data structure

	if (Accel_Struct == NONE) {  //no acceleration
		/* Get the closest intersection*/
		int num_objects = scene->getNumObjects();
		Object* obj = NULL;

		for (int i = 0; i < num_objects; i++)
		{
			obj = scene->getObject(i);
			auxRec = obj->hit(ray);    // Try intersecting the ray in an object

			/* If the ray intersects the nearest object then store both the object and its hit record */
			if (auxRec.isHit && auxRec.t < closestHit.t) {
				closestHit = auxRec;
				hitObj = obj;
			}
		}

		if (hitObj == NULL) {  // No intersected object
			if (skybox_flg)  //skybox cubemap overrides background color 
				color_Acc = scene->GetSkyboxColor(ray);
			else
				color_Acc = scene->GetBackgroundColor();

			return color_Acc.clamp();
		}
	}

	else if (Accel_Struct == GRID_ACC) {  // regular Grid
		if (!grid_ptr->Traverse(ray, &hitObj, closestHit)) {
			if (skybox_flg)
				color_Acc = scene->GetSkyboxColor(ray);
			else
				color_Acc = scene->GetBackgroundColor();
			return color_Acc.clamp();
		}
	}

	else if (Accel_Struct == BVH_ACC) { //BVH
		if (!bvh_ptr->Traverse(ray, &hitObj, closestHit)) {
			if (skybox_flg)
				color_Acc = scene->GetSkyboxColor(ray);
			else
				color_Acc = scene->GetBackgroundColor();
			return color_Acc.clamp();
		}
	}

	hitPoint = ray.origin + ray.direction * closestHit.t;
	N = closestHit.normal;

	// CALCULAR A COR DO PIXEL (Modelo de Reflexão Blinn-Phong)
	Material* mat = hitObj->GetMaterial(); // Obtém o material do objeto que foi atingido pelo raio
	Vector V = ray.direction * -1.0f;      // Vetor que aponta do ponto de interseção para a origem da câmera (visualização)

	// Percorre todas as luzes da cena para calcular a contribuição de cada uma
	for (int j = 0; j < num_lights; j++) {
		Light* light = scene->getLight(j); // Obtém a luz j da cena

		// Vetor do ponto atingido para a posição da luz (direção da luz)
		Vector L = (light->position - hitPoint).normalize();

		// Vetor "meio-termo" para o modelo Blinn-Phong: direção média entre a luz e a visão
		Vector H = (L + V).normalize();

		// ---------- Início Soft Shadows ----------
		 
		// Se AA ativado, usa spp, senão usa 4 amostras padrão
		const int numShadowSamples = (AA) ? spp : 4;

		Color diffuseTotal(0.0f, 0.0f, 0.0f);
		Color specularTotal(0.0f, 0.0f, 0.0f);

		// Loop para amostrar diferentes pontos na área da luz (para simular luzes com área e sombras suaves)
		for (int s = 0; s < numShadowSamples; s++) {
			Vector jitter; 

			if (AA) {
				// Se anti-aliasing está ativo, amostragem aleatória para suavizar bordas
				jitter = Vector(rand_float(), rand_float(), 0.0f);
			}
			else {
				// Se não, usa amostragem regular em grade (grid) 2x2 ou conforme numShadowSamples
				int i = s / numShadowSamples;
				int j = s % numShadowSamples;
				jitter = Vector((i + 0.5f) / numShadowSamples, (j + 0.5f) / numShadowSamples, 0.0f);
			}

			// Obtém um ponto amostrado na área da luz com base no jitter
			Vector lightPoint = light->getAreaLightPoint(jitter);

			// Vetor do ponto atingido para o ponto na área da luz (não normalizado ainda)
			Vector L = lightPoint - hitPoint;

			// Recalcula vetor "meio-termo" para este ponto específico da luz
			Vector H = (L + V).normalize();

			// Cria um raio sombra saindo do ponto de interseção na direção da luz
			Ray shadowRay(hitPoint + N * EPSILON, L);

			bool inShadow = false; 

			if (Accel_Struct == NONE) {
				// Se não tem estrutura de aceleração, faz teste de interseção com todos os objetos da cena
				for (int c = 0; c < scene->getNumObjects(); c++) {
					Object* shadowObj = scene->getObject(c);
					HitRecord shadowHit = shadowObj->hit(shadowRay);

					// Se houver interseção com objeto entre o ponto e a luz, ponto está em sombra
					if (shadowHit.isHit && shadowHit.t > EPSILON &&
						shadowHit.t < (lightPoint - hitPoint).length()) {
						inShadow = true;
						break;
					}
				}
			}
			else {
				if (Accel_Struct == GRID_ACC)
					inShadow = grid_ptr->Traverse(shadowRay);  
				else if (Accel_Struct == BVH_ACC)
					inShadow = bvh_ptr->Traverse(shadowRay);   
			}

			// Se ponto NÃO está em sombra, calcula contribuição difusa e especular para esta amostra
			if (!inShadow) {

				L = (lightPoint - hitPoint).normalize();

				// Produto escalar entre normal da superfície e direção da luz
				float NdotL = std::max(0.0f, N * L);

				// Produto escalar entre normal da superfície e vetor meio-termo para reflexão especular
				float NdotH = std::max(0.0f, N * H);

				// Contribuição difusa: cor difusa do material * luz * NdotL * coeficiente difuso do material
				diffuseTotal += mat->GetDiffColor() * light->emission * NdotL * mat->GetDiffuse();

				// Contribuição especular: cor especular do material * luz * (NdotH^brilho) * coeficiente especular
				specularTotal += mat->GetSpecColor() * light->emission * pow(NdotH, mat->GetShine()) * mat->GetSpecular();
			}
		}

		// Após somar todas as amostras, adiciona média das contribuições ao acumulador final de cor
		color_Acc += (diffuseTotal + specularTotal) * (1.0f / numShadowSamples);
	}

	// ------------------------------ Reflexão ------------------------------

	// Verifica se tem componente especular
	if (depth < MAX_DEPTH && mat->GetSpecular() > 0.0f) {
		// Calcula o vetor refletido da direção do raio em relação à normal N
		Vector reflectDir = ray.direction - N * 2.0f * (ray.direction * N);

		float roughness = 0.2f;

		// Gera um vetor aleatório dentro da esfera unitária para simular superfície rugosa
		Vector randomVec = random_in_unit_sphere();

		// Soma o vetor refletido com o vetor aleatório ponderado para dar um efeito "fuzzy"
		Vector fuzzyReflectDir = (reflectDir + randomVec * roughness).normalize();

		// Cria um novo raio refletido, com pequeno deslocamento na direção da normal para evitar acne de superfície
		Ray reflectRay(hitPoint + N * EPSILON, fuzzyReflectDir);

		// Chama recursivamente a função de ray-tracing para obter a cor da reflexão
		Color reflectColor = rayTracing(reflectRay, depth + 1, ior_1, lightSample);

		// Multiplica a cor refletida pela cor especular e coeficiente especular do material
		reflectColor *= mat->GetSpecColor() * mat->GetSpecular();

		// Adiciona a cor refletida ao acumulador final de cor do pixel
		color_Acc += reflectColor;
	}

	// ------------------------------ Refração com Fresnel ------------------------------

	// Se profundidade da recursão menor que limite e material tem componente de transmissão (transparência)
	if (depth < MAX_DEPTH && mat->GetTransmittance() > 0.0f) {
		Vector I = ray.direction; // Direção do raio de entrada
		Vector n = N;             // Normal da superfície no ponto de interseção

		// Índices de refração do meio de entrada e do material
		float eta_i = ior_1;
		float eta_t = mat->GetRefrIndex();

		// Calcula o cosseno do ângulo entre o raio e a normal, limitado entre -1 e 1
		float cos_i = std::clamp(I * N, -1.0f, 1.0f);

		// Verifica se o raio está dentro do objeto (inversão da normal e índices)
		bool insideObject = cos_i < 0.0f;

		if (insideObject) {
			// Se dentro do objeto, usa índices como estão e inverte o sinal do cosseno
			eta_i = ior_1;
			eta_t = mat->GetRefrIndex();
			n = N;
			cos_i = -cos_i;
		}
		else {
			// Se fora do objeto, inverte a normal e troca índices para simular saída do meio
			eta_i = ior_1;
			eta_t = 1.0f;
			n = N * -1.0f;
		}

		// Calcula razão dos índices de refração
		float eta = eta_i / eta_t;

		// Calcula parâmetro k para determinar se há refração real ou reflexão total interna
		float k = 1.0f - eta * eta * (1.0f - cos_i * cos_i);

		// Cálculo da reflexão usando aproximação de Schlick para o fator de Fresnel
		float R0 = pow((eta_i - eta_t) / (eta_i + eta_t), 2.0f);
		float k_r = R0 + (1.0f - R0) * pow(1.0f - cos_i, 5.0f);

		Color attenuation(1.0f, 1.0f, 1.0f); // Fator de atenuação (absorção)

		if (!insideObject) {
			// Se estamos fora do objeto, aplica a Lei de Beer para absorção da luz no material
			Color absorption = Color(1.0f, 1.0f, 1.0f) - mat->GetDiffColor();
			attenuation = (absorption * (-closestHit.t)).exp_(); // exp_() aplica exponencial para atenuação
		}
		else {
			// Se dentro do objeto, usa cor difusa do material sem atenuação
			attenuation = mat->GetDiffColor();
		}

		if (k >= 0.0f) {
			// Se k>=0, há refração (não houve reflexão total interna)
			Vector refractDir = (I * eta + n * (eta * cos_i - sqrtf(k))).normalize();

			// Cria raio refratado saindo do ponto de interseção com pequeno deslocamento para evitar acne
			Ray refractRay(hitPoint - n * EPSILON, refractDir);

			// Recursivamente calcula cor da refração, multiplicando por atenuação e (1 - fator de reflexão)
			Color refractColor = rayTracing(refractRay, depth + 1, eta_t, lightSample) * attenuation * (1.0f - k_r);

			color_Acc += refractColor;

		}
		else {
			// Caso contrário, ocorre reflexão total interna, o fator de reflexão é 1
			k_r = 1.0f;
		}

		// Cálculo da contribuição da reflexão usando o fator de Fresnel k_r
		Vector reflectDir = I - n * 2.0f * (I * n);
		reflectDir.normalize();

		// Cria raio refletido para o lado oposto com pequeno deslocamento para evitar acne
		Ray reflectRay(hitPoint + n * EPSILON, reflectDir);

		// Calcula recursivamente a cor da reflexão, ponderada pelo fator k_r e cor especular do material
		Color reflectColor = rayTracing(reflectRay, depth + 1, eta_i, lightSample) * k_r * mat->GetSpecColor();

		color_Acc += reflectColor;
	}

	// Retorna a cor final do pixel, limitando os valores para o intervalo válido (0 a 1)
	return color_Acc.clamp();
}


// Render function by primary ray casting from the eye towards the scene's objects
void renderScene()	
{
	unsigned int counter = 0;
	set_rand_seed(time(NULL) * time(NULL)); // Use current time as seed for random generator 

	if (drawModeEnabled) {
		//glClear(GL_COLOR_BUFFER_BIT);
		scene->GetCamera()->SetEye(Vector(camX, camY, camZ));  //Camera motion
	}
	

	if (Progressive_flg){          ///////////////////////// ZONE A  - Progressive RayTracer/////////////////////////////////////////
		if (FrameCount < MAX_SAMPLES) {
#pragma omp parallel for collapse(2)
			for (int y = 0; y < RES_Y; y++) {
				for (int x = 0; x < RES_X; x++) {

					int index_pos = 0;
					int index_col = 0;
					Color color;
					Ray ray;
					Vector pixel_sample;  //viewport coordinates
					Vector light_sample = Vector(0.0f, 0.0f, 0.0f); // sample in Light coordinates

					pixel_sample.x = x + rand_double();
					pixel_sample.y = y + rand_double();

					if (!DOF) ray = scene->GetCamera()->PrimaryRay(pixel_sample);
					else {        // sample_unit_disk() returns [-1 1] and aperture is the diameter of the lens

						Vector lens_sample = rnd_unit_disk() * scene->GetCamera()->GetAperture() / 2.0f;  // lens sample in Camera coordinates

						ray = scene->GetCamera()->PrimaryRay(lens_sample, pixel_sample);
					}
					color = rayTracing(ray, 1, 1.0, Vector(rand_float(), rand_float(), 0.0f));

					index_pos = 2 * (x + RES_X * y);
					vertices[index_pos] = (float)x;
					vertices[index_pos + 1] = (float)y;

					index_col = 3 * (x + RES_X * y);
					///////////UNDERSTAND AND EXPLAIN THE FOLLOWING CODE//////////////////////
					if (FrameCount == 1) {
						colors[index_col] = (float)color.r();
						colors[index_col + 1] = (float)color.g();
						colors[index_col + 2] = (float)color.b();
					}

					else {
						colors[index_col] = lerp(colors[index_col], (float)color.r(), 1.0 / FrameCount);
						index_col++;
						colors[index_col] = lerp(colors[index_col], (float)color.g(), 1.0 / FrameCount);
						index_col++;
						colors[index_col] = lerp(colors[index_col], (float)color.b(), 1.0 / FrameCount);
					}
				}
			}
		}
		drawPoints();
		if(FrameCount != MAX_SAMPLES)  FrameCount++;
		FPS++;  
		std::ostringstream oss;
		oss << CAPTION << ": " << FramesPerSecond << " FPS @ (" << RES_X << "x" << RES_Y << ") @ Samples number: " << FrameCount;
		std::string s = oss.str();
		glutSetWindow(WindowHandle);
		glutSetWindowTitle(s.c_str());
		glutSwapBuffers();
	}

	//////////////ZONE B - NOT Progressive RayTracer////////////////////////////////////
	else {				
#pragma omp parallel for collapse(2)
		for (int y = 0; y < RES_Y; y++) {
			for (int x = 0; x < RES_X; x++) {
				Color color;
				int index_pos = 0;
				int index_col = 0;
				Ray ray;
				Vector pixel_sample;  //viewport coordinates
				Vector light_sample; // sample in Light coordinates

				// ZONE B.1 - Distribution Ray Tracer: pixel, area light and lens supersampling with jittering
				if (AA) {
					Color temp_color(0.0f, 0.0f, 0.0f);
					float r_sum = 0.0f, g_sum = 0.0f, b_sum = 0.0f;

					#pragma omp parallel for reduction(+:r_sum, g_sum, b_sum)
					for (int p = 0; p < spp; p++) {
						float jitter_x = rand_float();
						float jitter_y = rand_float();
						Vector pixel_sample(x + jitter_x, y + jitter_y, 0.0f);

						Vector lens_sample(0.0f);
						if (DOF) {
							lens_sample = rnd_unit_disk() * scene->GetCamera()->GetAperture() / 2.0f;
						}

						Ray ray = DOF ? scene->GetCamera()->PrimaryRay(lens_sample, pixel_sample)
							: scene->GetCamera()->PrimaryRay(pixel_sample);

						// Light jittering (soft shadow)
						Vector light_sample = Vector(rand_float(), rand_float(), 0.0f);

						Color sample = rayTracing(ray, 1, 1.0, light_sample);
						r_sum += sample.r();
						g_sum += sample.g();
						b_sum += sample.b();
					}

					temp_color = Color(r_sum, g_sum, b_sum);
					color = temp_color * (1.0f / spp);
				}

				//ZONE B.2 - Whitted ray tracer (without antialiasing)
				else {	
					pixel_sample.x = x + 0.5f;  
					pixel_sample.y = y + 0.5f;
					Ray ray1 = scene->GetCamera()->PrimaryRay(pixel_sample);

					Color accumulated_color(0.0f, 0.0f, 0.0f);

					// Number of samples per area light (n x n)
					const int n = 2;

					for (int i = 0; i < n; ++i) {
						for (int j = 0; j < n; ++j) {
							Vector light_sample((i + 0.5f) / n, (j + 0.5f) / n, 0.0f);
							accumulated_color += rayTracing(ray1, 1, 1.0, light_sample);
						}
					}

					color = accumulated_color * (1.0f / (n * n));
				}

				if (drawModeEnabled) {
					index_pos = 2 * (x + RES_X * y);
					vertices[index_pos] = (float)x;
					vertices[index_pos + 1] = (float)y;

					index_col = 3 * (x + RES_X * y);
					colors[index_col] = (float)color.r();
					colors[index_col+1] = (float)color.g();
					colors[index_col+2] = (float)color.b();
				}
				else {
					img_Data[counter++] = u8fromfloat((float)color.r());
					img_Data[counter++] = u8fromfloat((float)color.g());
					img_Data[counter++] = u8fromfloat((float)color.b());
				}
			}
		}

		if (drawModeEnabled){
			FPS++;
			drawPoints();
			glutSwapBuffers();
		}
		else {
			printf("Rendering ended!\n"); 	
			if (saveImgFile("RT_Output.png") != IL_NO_ERROR) {
				printf("Error saving Image file\n");
				exit(0);
			}
			printf("Image file created\n");
		}
	}
}

// Callback function for glutCloseFunc
void cleanup()
{
	destroyShaderProgram();
	destroyBufferObjects();
}

void ortho(float left, float right, float bottom, float top, 
			float nearp, float farp)
{
	m[0 * 4 + 0] = 2 / (right - left);
	m[0 * 4 + 1] = 0.0;
	m[0 * 4 + 2] = 0.0;
	m[0 * 4 + 3] = 0.0;
	m[1 * 4 + 0] = 0.0;
	m[1 * 4 + 1] = 2 / (top - bottom);
	m[1 * 4 + 2] = 0.0;
	m[1 * 4 + 3] = 0.0;
	m[2 * 4 + 0] = 0.0;
	m[2 * 4 + 1] = 0.0;
	m[2 * 4 + 2] = -2 / (farp - nearp);
	m[2 * 4 + 3] = 0.0;
	m[3 * 4 + 0] = -(right + left) / (right - left);
	m[3 * 4 + 1] = -(top + bottom) / (top - bottom);
	m[3 * 4 + 2] = -(farp + nearp) / (farp - nearp);
	m[3 * 4 + 3] = 1.0;
}

void reshape(int w, int h)
{
    glClear(GL_COLOR_BUFFER_BIT);
	glViewport(0, 0, w, h);
	ortho(0, (float)RES_X, 0, (float)RES_Y, -1.0, 1.0);
}

void processKeys(unsigned char key, int xx, int yy)
{
	switch (key) {

	case 27:
		glutLeaveMainLoop();
		break;

	case 'p':
		if (Progressive_flg) Progressive_flg = false;
		else { Progressive_flg = true; FrameCount = 1; }
		break;

	case 'r':
		camX = Eye.x;
		camY = Eye.y;
		camZ = Eye.z;
		r = Eye.length();
		_beta = asinf(camY / r) * 180.0f / 3.14f;
		alpha = atanf(camX / camZ) * 180.0f / 3.14f;
		FrameCount = 1;
		break;

	case 'c':
		printf("Camera Spherical Coordinates (%f, %f, %f)\n", r, _beta, alpha);
		printf("Camera Cartesian Coordinates (%f, %f, %f)\n", camX, camY, camZ);
		break;
	}
}


// ------------------------------------------------------------
//
// Mouse Events
//

void processMouseButtons(int button, int state, int xx, int yy)
{
	// start tracking the mouse
	if (state == GLUT_DOWN) {
		startX = xx;
		startY = yy;
		//FrameCount = 1;
		if (button == GLUT_LEFT_BUTTON)
			tracking = 1;
		else if (button == GLUT_RIGHT_BUTTON)
			tracking = 2;
	}

	//stop tracking the mouse
	else if (state == GLUT_UP) {
		if (tracking == 1) {
			alpha -= (xx - startX);
			_beta += (yy - startY);
		}
		else if (tracking == 2) {
			r += (yy - startY) * 0.01f;
			if (r < 0.1f)
				r = 0.1f;
		}
		tracking = 0;
	}
}

// Track mouse motion while buttons are pressed

void processMouseMotion(int xx, int yy)
{

	int deltaX, deltaY;
	float alphaAux, betaAux;
	float rAux;

	deltaX = -xx + startX;
	deltaY = yy - startY;

	// left mouse button: move camera
	if (tracking == 1) {


		alphaAux = alpha + deltaX;
		betaAux = _beta + deltaY;

		if (betaAux > 85.0f)
			betaAux = 85.0f;
		else if (betaAux < -85.0f)
			betaAux = -85.0f;
		rAux = r;
	}
	// right mouse button: zoom
	else if (tracking == 2) {

		alphaAux = alpha;
		betaAux = _beta;
		rAux = r + (deltaY * 0.01f);
		if (rAux < 0.1f)
			rAux = 0.1f;
	}

	FrameCount = 1;
	camX = rAux * sin(alphaAux * 3.14f / 180.0f) * cos(betaAux * 3.14f / 180.0f);
	camZ = rAux * cos(alphaAux * 3.14f / 180.0f) * cos(betaAux * 3.14f / 180.0f);
	camY = rAux * sin(betaAux * 3.14f / 180.0f);

	//  uncomment this if not using an idle func
	//	glutPostRedisplay();
}

void mouseWheel(int wheel, int direction, int x, int y) {

	r += direction * 0.1f;
	if (r < 0.1f)
		r = 0.1f;
	FrameCount = 1;
	camX = r * sin(alpha * 3.14f / 180.0f) * cos(_beta * 3.14f / 180.0f);
	camZ = r * cos(alpha * 3.14f / 180.0f) * cos(_beta * 3.14f / 180.0f);
	camY = r * sin(_beta * 3.14f / 180.0f);

	//  uncomment this if not using an idle func
	//	glutPostRedisplay();
}


/////////////////////////////////////////////////////////////////////// SETUP

void setupCallbacks() 
{
	glutKeyboardFunc(processKeys);
	glutCloseFunc(cleanup);
	glutDisplayFunc(renderScene);
	glutReshapeFunc(reshape);
	glutMouseFunc(processMouseButtons);
	glutMotionFunc(processMouseMotion);
	glutMouseWheelFunc(mouseWheel);
	glutTimerFunc(1000, timer, 0);
	glutIdleFunc(renderScene);
}

void setupGLEW() {
	glewExperimental = GL_TRUE;
	GLenum result = glewInit() ; 
	if (result != GLEW_OK) { 
		std::cerr << "ERROR glewInit: " << glewGetString(result) << std::endl;
		exit(EXIT_FAILURE);
	} 
	GLenum err_code = glGetError();
	printf ("Vendor: %s\n", glGetString (GL_VENDOR));
	printf ("Renderer: %s\n", glGetString (GL_RENDERER));
	printf ("Version: %s\n", glGetString (GL_VERSION));
	printf ("GLSL: %s\n", glGetString (GL_SHADING_LANGUAGE_VERSION));
}

void setupGLUT(int argc, char* argv[])
{
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA);
	glutInitContextVersion(4, 3);
	glutInitContextFlags(GLUT_FORWARD_COMPATIBLE);
	glutInitContextProfile(GLUT_CORE_PROFILE);
	glutSetOption(GLUT_ACTION_ON_WINDOW_CLOSE,GLUT_ACTION_GLUTMAINLOOP_RETURNS);
	glutInitWindowPosition(100,250);
	glutInitWindowSize(RES_X, RES_Y);
	
	glDisable(GL_DEPTH_TEST);
	WindowHandle = glutCreateWindow(CAPTION);
	if(WindowHandle < 1) {
		std::cerr << "ERROR: Could not create a new rendering window." << std::endl;
		exit(EXIT_FAILURE);
	}
}


void init(int argc, char* argv[])
{
	// set the initial camera position on its spherical coordinates
	Eye  =  scene->GetCamera()->GetEye();
	camX = Eye.x;
	camY = Eye.y;
	camZ = Eye.z;
	r = Eye.length();
	_beta = asinf(camY/r) * 180.0f / 3.14f;
	alpha = atanf(camX / camZ) * 180.0f / 3.14f;

	setupGLUT(argc, argv);
	setupCallbacks();
	setupGLEW();
	glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
	std::cerr << "CONTEXT: OpenGL v" << glGetString(GL_VERSION) << std::endl;
	createShaderProgram();
	createBufferObjects();
}

void init_scene(void)
{
	char scenes_dir[70] = "P3D_Scenes/";
	char input_user[50];
	char scene_name[70];

	scene = new Scene();

	if (P3F_scene) {  //Loading a P3F scene

		while (true) {
			cout << "Input the Scene Name: ";
			cin >> input_user;
			strcpy_s(scene_name, sizeof(scene_name), scenes_dir);
			strcat_s(scene_name, sizeof(scene_name), input_user);

			ifstream file(scene_name, ios::in);
			if (file.fail()) {
				printf("\nError opening P3F file.\n");
			}
			else
				break;
		}

		scene->load_p3f(scene_name);
		printf("Scene loaded.\n\n");
	}
	else {
		printf("Creating Peter Shirley Scene.\n\n");
		scene->create_random_scene();
	}

	RES_X = scene->GetCamera()->GetResX();
	RES_Y = scene->GetCamera()->GetResY();
	printf("\nResolutionX = %d  ResolutionY= %d.\n", RES_X, RES_Y);

	spp = scene->GetSamplesPerPixel();
	if (spp != 0) {
		AA = true; //Anti-aliasing
		printf("\nJittering with %d samples per pixel.\n", scene->GetSamplesPerPixel());
	}
	else AA = false;

	
	if ((scene->GetCamera()->GetAperture() != 0) && AA) {
		DOF = true; //Depth-Of-Field effect enabled only if AA enabled
		printf("Depth-Of-Field effect enabled\n");
	}
	else DOF = false;

	// Pixel buffer to be used in the Save Image function
	img_Data = (uint8_t*)malloc(3 * RES_X*RES_Y * sizeof(uint8_t));
	if (img_Data == NULL) exit(1);

	Accel_Struct = scene->GetAccelStruct();   //Type of acceleration data structure

	if (Accel_Struct == GRID_ACC) {
		grid_ptr = new Grid();
		vector<Object*> objs;
		int num_objects = scene->getNumObjects();

		for (int o = 0; o < num_objects; o++) {
			objs.push_back(scene->getObject(o));
		}
		grid_ptr->Build(objs);
		printf("Grid built.\n\n");
	}
	else if (Accel_Struct == BVH_ACC) {
		vector<Object*> objs;
		int num_objects = scene->getNumObjects();
		bvh_ptr = new BVH();

		for (int o = 0; o < num_objects; o++) {
			objs.push_back(scene->getObject(o));
		}
		bvh_ptr->Build(objs);
		printf("BVH built.\n\n");
	}
	else
		printf("No acceleration data structure.\n\n");

	unsigned int spp = scene->GetSamplesPerPixel();
	if (spp == 0)
		printf("Whitted Ray-Tracing\n");
	else
		printf("Distribution Ray-Tracing\n");
}

int main(int argc, char* argv[])
{
	//Initialization of DevIL 
	if (ilGetInteger(IL_VERSION_NUM) < IL_VERSION)
	{
		printf("wrong DevIL version \n");
		exit(0);
	}
	ilInit();

	int ch;
	if (!drawModeEnabled) {

		do {
			init_scene();
			
			auto timeStart = std::chrono::high_resolution_clock::now();
			renderScene();  //Just creating an image file
			auto timeEnd = std::chrono::high_resolution_clock::now();
			auto passedTime = std::chrono::duration<double, std::milli>(timeEnd - timeStart).count();
			printf("\nDone: %.2f (sec)\n", passedTime / 1000);
			if (!P3F_scene) break;
			cout << "\nPress 'y' to render another image or another key to terminate!\n";
			delete(scene);
			if (Accel_Struct == GRID_ACC) delete(grid_ptr);
			else if (Accel_Struct == BVH_ACC) delete(bvh_ptr);
			free(img_Data);
			ch = _getch();
		} while((toupper(ch) == 'Y')) ;
	}

	else {   //Use OpenGL to draw image in the screen
		printf("OPENGL DRAWING MODE\n\n");
		init_scene();
		
	//	if (Accel_Struct == 1) grid_ptr->Build();  
		
		size_vertices = 2 * RES_X * RES_Y * sizeof(float);
		size_colors = 3 * RES_X * RES_Y * sizeof(float);
		vertices = (float*)malloc(size_vertices);
		if (vertices == NULL) exit(1);
		colors = (float*)malloc(size_colors);
		if (colors == NULL) exit(1);
		memset(colors, 0, size_colors);
		/* Setup GLUT and GLEW */
		init(argc, argv);
		glutMainLoop();
	}

	free(colors);
	free(vertices);
	printf("Programa terminado normalmente\n");
	exit(EXIT_SUCCESS);
}
///////////////////////////////////////////////////////////////////////