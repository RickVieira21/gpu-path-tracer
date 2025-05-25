#include "rayAccelerator.h"
#include "macros.h"

using namespace std;

BVH::BVHNode::BVHNode(void) {}

void BVH::BVHNode::setAABB(AABB& bbox_) { this->bbox = bbox_; }

void BVH::BVHNode::makeLeaf(unsigned int index_, unsigned int n_objs_) {
	this->leaf = true;
	this->index = index_; 
	this->n_objs = n_objs_; 
}

void BVH::BVHNode::makeNode(unsigned int left_index_) {
	this->leaf = false;
	this->index = left_index_; 
}


BVH::BVH(void) {}

int BVH::getNumObjects() { return objects.size(); }


void BVH::Build(vector<Object *> &objs) {

		
			BVHNode *root = new BVHNode();

			Vector min = Vector(FLT_MAX, FLT_MAX, FLT_MAX), max = Vector(-FLT_MAX, -FLT_MAX, -FLT_MAX);
			AABB world_bbox = AABB(min, max);

			for (Object* obj : objs) {
				AABB bbox = obj->GetBoundingBox();
				world_bbox.extend(bbox);
				objects.push_back(obj);
			}
			world_bbox.min.x -= EPSILON; world_bbox.min.y -= EPSILON; world_bbox.min.z -= EPSILON;
			world_bbox.max.x += EPSILON; world_bbox.max.y += EPSILON; world_bbox.max.z += EPSILON;
			root->setAABB(world_bbox);
			nodes.push_back(root);
			build_recursive(0, objects.size(), root); // -> root node takes all the objects
		}





// BVH::build_recursive - Constrói recursivamente a BVH 
// 1 - A função divide recursivamente os objetos num nó interior com dois filhos (esquerdo e direito),
// com base no centroide dos objetos ao longo do eixo mais longo do bounding box atual.
// 2 - Quando o número de objetos é menor ou igual a um limiar (threshold), cria um nó folha.
// 3 - Os nós filhos são adicionados à lista global de nós (`nodes`) e os seus bounding boxes são computados.
// 4 - A estrutura final é uma árvore binária onde cada nó representa um AABB (Axis-Aligned Bounding Box) que engloba
// os seus objetos ou filhos.

void BVH::build_recursive(int left_index, int right_index, BVHNode* node) {

	const int threshold = 2; // Limite mínimo de objetos por folha. 

	int n_objects = right_index - left_index; // Número de objetos na sub-região atual do vetor de objetos

	// Caso base: se o número de objetos for igual ou inferior ao limiar, criamos um nó folha.
	if (n_objects <= threshold) {
		node->makeLeaf(left_index, n_objects); // Configura o nó como folha e associa os objetos.
		return;
	}

	// Calcular a Bounding Box (AABB) que engloba todos os objetos entre left_index e right_index.
	AABB bbox;
	for (int i = left_index; i < right_index; ++i) {
		bbox.extend(objects[i]->GetBoundingBox()); // Expande a caixa com a AABB do objeto atual
	}

	// Determinar o eixo mais longo da caixa (x = 0, y = 1, z = 2)
	// Este será o eixo onde faremos a divisão dos objetos (maior dispersão)
	Vector diag = bbox.max - bbox.min;
	int axis = 0;
	if (diag.y > diag.x && diag.y > diag.z) axis = 1;
	else if (diag.z > diag.x && diag.z > diag.y) axis = 2;

	// Ordenar os objetos de acordo com o centroide no eixo escolhido
	auto centroidComparator = [axis](Object* a, Object* b) {
		const Vector centroidA = a->getCentroid();
		const Vector centroidB = b->getCentroid();
		float valA, valB;

		switch (axis) {
		case 0: valA = centroidA.x; valB = centroidB.x; break;
		case 1: valA = centroidA.y; valB = centroidB.y; break;
		default: valA = centroidA.z; valB = centroidB.z; break;
		}

		return valA < valB; // Ordena de forma crescente
		};

	std::sort(objects.begin() + left_index, objects.begin() + right_index, centroidComparator);

	// Encontrar o ponto de divisão (split) baseado na média dos centroides no eixo escolhido
	auto getCoord = [](const Vector& v, int axis) -> float {
		switch (axis) {
		case 0: return v.x;
		case 1: return v.y;
		case 2: return v.z;
		default: return 0.0f;
		}
		};

	float split_coord = 0.0f;
	for (int i = left_index; i < right_index; ++i) {
		split_coord += getCoord(objects[i]->GetBoundingBox().centroid(), axis);
	}
	split_coord /= n_objects; // Média

	// Determina split_index: o primeiro índice cujo centroide é maior que a média
	int split_index = left_index;
	for (; split_index < right_index; ++split_index) {
		if (getCoord(objects[split_index]->GetBoundingBox().centroid(), axis) > split_coord)
			break;
	}

	// Garante que nenhuma partição fique vazia (caso todos objetos fiquem de um lado só)
	if (split_index == left_index || split_index == right_index) {
		split_index = left_index + (n_objects / 2); // Divisão ao meio
	}

	// Criar os nós filho (esquerdo e direito)
	BVHNode* leftNode = new BVHNode();
	BVHNode* rightNode = new BVHNode();

	// Calcular AABB para os filhos
	AABB left_bbox, right_bbox;
	for (int i = left_index; i < split_index; ++i)
		left_bbox.extend(objects[i]->GetBoundingBox());
	for (int i = split_index; i < right_index; ++i)
		right_bbox.extend(objects[i]->GetBoundingBox());

	leftNode->setAABB(left_bbox);
	rightNode->setAABB(right_bbox);

	// Armazenar os nós filhos no vetor global `nodes`
	int left_node_index = nodes.size();
	nodes.push_back(leftNode);   // Esquerdo
	nodes.push_back(rightNode);  // Direito

	// Configura o nó atual como nó interno apontando para o filho esquerdo
	node->makeNode(left_node_index);

	// Chamada recursiva para construir os filhos
	build_recursive(left_index, split_index, leftNode);       // esquerda
	build_recursive(split_index, right_index, rightNode);     // direita
}







// BVH::Traverse (para objetos visiveis) - Percorre a árvore BVH com um raio de interseção completo, procurando o objeto mais próximo que colide com o raio.
// 1 - Utiliza uma `hit_stack` para evitar recursão, navegando nos nós da BVH do mais próximo para o mais distante, com base na interseção do raio com os bounding boxes dos nós.
// 2 - Quando encontra um nó folha, verifica interseção com todos os objetos contidos nele.
// 3 - Se encontrar uma interseção válida e mais próxima que as anteriores, atualiza o objeto atingido `hit_obj` e o `hitRec` correspondente.
// Retorna `true` se algum objeto for atingido pelo raio; caso contrário, retorna `false`.

bool BVH::Traverse(Ray& ray, Object** hit_obj, HitRecord& hitRec) {
	float t_closest = FLT_MAX;          // Guarda o t mais próximo encontrado
	bool hit = false;                   
	stack<StackItem> hit_stack;         
	BVHNode* currentNode = nodes[0];    // Começa pela raiz da árvore

	float root_t = FLT_MAX;

	// Se o raio não intersecta a AABB da raiz, não há colisão
	if (!currentNode->getAABB().hit(ray, root_t))
		return false;

	while (true) {
		if (!currentNode->isLeaf()) {
			// Nó interno: obter filhos esquerdo e direito
			BVHNode* left = nodes[currentNode->getIndex()];
			BVHNode* right = nodes[currentNode->getIndex() + 1];

			float t_left = FLT_MAX, t_right = FLT_MAX;
			bool hit_left = left->getAABB().hit(ray, t_left);
			bool hit_right = right->getAABB().hit(ray, t_right);

			if (hit_left && hit_right) {
				// Ambos filhos atingidos: prioriza o mais próximo, empilha o outro
				if (t_left < t_right) {
					hit_stack.push(StackItem(right, t_right));
					currentNode = left;
				}
				else {
					hit_stack.push(StackItem(left, t_left));
					currentNode = right;
				}
			}
			else if (hit_left) currentNode = left;
			else if (hit_right) currentNode = right;
			else goto pop_stack;
		}
		else {
			// Nó folha: verificar colisões com cada objeto
			unsigned int start = currentNode->getIndex();
			unsigned int count = currentNode->getNObjs();

			for (unsigned int i = 0; i < count; ++i) {
				HitRecord tempRec = objects[start + i]->hit(ray);
				if (tempRec.isHit && tempRec.t < t_closest) {
					t_closest = tempRec.t;
					hitRec = tempRec;
					*hit_obj = objects[start + i];
					hit = true;
				}
			}

		pop_stack:
			// Retomar nós empilhados
			while (!hit_stack.empty()) {
				StackItem item = hit_stack.top();
				hit_stack.pop();
				if (item.t < t_closest) {
					currentNode = item.ptr;
					goto continue_loop;
				}
			}
			break; // stack vazia: fim do traverse
		}

	continue_loop:;
	}

	return hit; // true se algum objeto foi atingido
}






// BVH::Traverse (shadow rays) - Percorre a árvore BVH com um shadow ray (raio de sombra) para determinar se há obstruções entre um ponto e a luz.
// 1 - Normaliza a direção do raio e limita a sua extensão até à distância da luz (ray length).
// 2 - Navega pela BVH de forma semelhante ao Traverse principal, mas retorna `true` assim que encontrar qualquer interseção válida
// Retorna `true` se houver algum objeto entre o ponto e a luz (há sombra), e `false` se o caminho estiver livre.

bool BVH::Traverse(Ray& ray) {
	stack<StackItem> hit_stack;

	double length = ray.direction.length(); // Distância até a fonte de luz
	Vector dir = ray.direction.normalize();
	Ray shadowRay(ray.origin, dir); // Ray normalizado para uso local

	BVHNode* currentNode = nodes[0];
	float root_t = FLT_MAX;

	if (!currentNode->getAABB().hit(shadowRay, root_t))
		return false; // Nenhuma interseção com a raiz

	while (true) {
		if (!currentNode->isLeaf()) {
			BVHNode* left = nodes[currentNode->getIndex()];
			BVHNode* right = nodes[currentNode->getIndex() + 1];

			float t_left = FLT_MAX, t_right = FLT_MAX;
			bool hit_left = left->getAABB().hit(shadowRay, t_left);
			bool hit_right = right->getAABB().hit(shadowRay, t_right);

			if (hit_left && hit_right) {
				if (t_left < t_right) {
					hit_stack.push(StackItem(right, t_right));
					currentNode = left;
				}
				else {
					hit_stack.push(StackItem(left, t_left));
					currentNode = right;
				}
			}
			else if (hit_left) currentNode = left;
			else if (hit_right) currentNode = right;
			else goto pop_stack;
		}
		else {
			unsigned int start = currentNode->getIndex();
			unsigned int count = currentNode->getNObjs();

			for (unsigned int i = 0; i < count; ++i) {
				HitRecord tempRec = objects[start + i]->hit(shadowRay);

				// Qualquer colisão válida dentro da distância para a luz bloqueia a luz (sombra)
				if (tempRec.isHit && tempRec.t > EPSILON && tempRec.t < length)
					return true; // Está na sombra
			}

		pop_stack:
			while (!hit_stack.empty()) {
				StackItem item = hit_stack.top();
				hit_stack.pop();
				if (item.t < length) {
					currentNode = item.ptr;
					goto continue_loop;
				}
			}

			break; // stack vazia, fim
		}

	continue_loop:;
	}

	return false; // Sem interseção: está iluminado
}

