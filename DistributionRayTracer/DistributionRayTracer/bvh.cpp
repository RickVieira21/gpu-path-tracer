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

void BVH::build_recursive(int left_index, int right_index, BVHNode *node) {

	//PUT YOUR CODE HERE

		//right_index, left_index and split_index refer to the indices in the objects vector
	   // do not confuse with left_nodde_index and right_node_index which refer to indices in the nodes vector. 
	    // node.index can have a index of objects vector or a index of nodes vector

	const int threshold = 2; // Pode ajustar dependendo da aplicação

	int n_objects = right_index - left_index;

	// Caso base: criar folha
	if (n_objects <= threshold) {
		node->makeLeaf(left_index, n_objects);
		return;
	}

	// Calcular bounding box da região atual
	AABB bbox;
	for (int i = left_index; i < right_index; ++i) {
		bbox.extend(objects[i]->GetBoundingBox());
	}

	// Determinar o maior eixo
	Vector diag = bbox.max - bbox.min;
	int axis = 0;
	if (diag.y > diag.x && diag.y > diag.z) axis = 1;
	else if (diag.z > diag.x && diag.z > diag.y) axis = 2;

	// Ordenar objetos de acordo com o centro no maior eixo
	auto centroidComparator = [axis](Object* a, Object* b) {
		const Vector centroidA = a->getCentroid();
		const Vector centroidB = b->getCentroid();
		float valA, valB;

		switch (axis) {
		case 0:
			valA = centroidA.x;
			valB = centroidB.x;
			break;
		case 1:
			valA = centroidA.y;
			valB = centroidB.y;
			break;
		case 2:
		default:
			valA = centroidA.z;
			valB = centroidB.z;
			break;
		}

		return valA < valB;
		};

	std::sort(objects.begin() + left_index, objects.begin() + right_index, centroidComparator);

	// Encontrar índice de divisão pelo centroide médio
	auto getCoord = [](const Vector& v, int axis) -> float {
		switch (axis) {
			case 0: return v.x;
			case 1: return v.y;
			case 2: return v.z;
			default: return 0.0f; // fallback de segurança
		}
	};

	float split_coord = 0.0f;
	for (int i = left_index; i < right_index; ++i) {
		split_coord += getCoord(objects[i]->GetBoundingBox().centroid(), axis);
	}

	split_coord /= n_objects;
	int split_index = left_index;

	for (; split_index < right_index; ++split_index) {
		if (getCoord(objects[split_index]->GetBoundingBox().centroid(), axis) > split_coord) break;
	}

	// Evitar partições vazias
	if (split_index == left_index || split_index == right_index) {
		split_index = left_index + (n_objects / 2);
	}

	// Criar nós esquerdo e direito
	BVHNode* leftNode = new BVHNode();
	BVHNode* rightNode = new BVHNode();

	AABB left_bbox, right_bbox;
	for (int i = left_index; i < split_index; ++i)
		left_bbox.extend(objects[i]->GetBoundingBox());
	for (int i = split_index; i < right_index; ++i)
		right_bbox.extend(objects[i]->GetBoundingBox());

	leftNode->setAABB(left_bbox);
	rightNode->setAABB(right_bbox);

	// Inserir filhos no vetor de nós e configurar este nó como interior
	int left_node_index = nodes.size();
	nodes.push_back(leftNode);
	nodes.push_back(rightNode);
	node->makeNode(left_node_index); // index aponta para o filho esquerdo

	// Chamadas recursivas
	build_recursive(left_index, split_index, leftNode);
	build_recursive(split_index, right_index, rightNode);
			
	}

bool BVH::Traverse(Ray& ray, Object** hit_obj, HitRecord& hitRec) {
			float tmp;
			bool hit = false;
			stack<StackItem> hit_stack;
			HitRecord rec;   //rec.isHit initialized to false and rec.t initialized with FLT_MAX

			BVHNode* currentNode = nodes[0];

			//PUT YOUR CODE HERE
			float t_closest = FLT_MAX;
			float root_t = FLT_MAX;

			if (!currentNode->getAABB().hit(ray, root_t)) {
				return false; // Raio não atinge a AABB da raiz
			}

			while (true) {
				if (!currentNode->isLeaf()) {
					BVHNode* left = nodes[currentNode->getIndex()];
					BVHNode* right = nodes[currentNode->getIndex() + 1];

					float t_left = FLT_MAX, t_right = FLT_MAX;
					bool hit_left = left->getAABB().hit(ray, t_left);
					bool hit_right = right->getAABB().hit(ray, t_right);

					if (hit_left && hit_right) {
						// Ambos atingidos: empilha o mais distante
						if (t_left < t_right) {
							hit_stack.push(StackItem(right, t_right));
							currentNode = left;
						}
						else {
							hit_stack.push(StackItem(left, t_left));
							currentNode = right;
						}
					}
					else if (hit_left) {
						currentNode = left;
					}
					else if (hit_right) {
						currentNode = right;
					}
					else {
						// Nenhuma interseção: prossegue na pilha
						goto pop_stack;
					}
				}
				else {
					// Nó folha: testar interseções com os objetos
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

					// Avança na pilha
				pop_stack:
					while (!hit_stack.empty()) {
						StackItem item = hit_stack.top();
						hit_stack.pop();
						if (item.t < t_closest) {
							currentNode = item.ptr;
							goto continue_loop;
						}
					}
					break; // Pilha vazia, finaliza
				}

			continue_loop:;
			}

			return hit;
		
	}

bool BVH::Traverse(Ray& ray) {  //shadow ray with length
			float tmp;
			stack<StackItem> hit_stack;
			HitRecord rec;

			double length = ray.direction.length(); //distance between light and intersection point
			ray.direction.normalize();

		
			return false;  //no primitive intersection
			
	}		
