#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <iostream>
#include <fstream>
#include <cmath>
#include <ctime>
#include <cstdlib>
#include <map>
#include <vector>
#include <queue>
typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef CGAL::Polyhedron_3<Kernel> Polyhedron;
typedef CGAL::Vector_3<Kernel> Vector3;
typedef Kernel::Point_3 Point_3;
typedef Polyhedron::Facet_iterator Facet_iterator;
typedef Polyhedron::Vertex_iterator Vertex_iterator;
typedef Polyhedron::Halfedge_iterator Halfedge_iterator;
typedef Polyhedron::Halfedge_around_facet_circulator Halfedge_facet_circulator;
typedef std::map<Polyhedron::Facet_handle, double> Facet_double_map;
typedef std::map<Polyhedron::Facet_handle, int> Facet_int_map;
typedef std::map<Polyhedron::Facet_handle, bool> Facet_bool_map;
typedef std::vector<double> Vect_double;
typedef std::vector<int> Vect_int;
typedef std::vector<Vect_double> Vect_Color;
typedef std::queue<Polyhedron::Facet_handle> Facet_Queue;

double maxMeasure = 0.0;
double minMeasure = DBL_MAX;
double moyMeasure = 0.0;

double maxArea = 0.0;
double minArea = DBL_MAX;

void calculPerimetre(Polyhedron &mesh, Facet_double_map &map) {
	Facet_iterator f_it = mesh.facets_begin();
	int nbFacets = 0;
	while (f_it != mesh.facets_end()) {
		Halfedge_facet_circulator h_it = f_it->facet_begin();
		CGAL_assertion(CGAL::circulator_size(h_it) >= 3);
		double perim = 0.0;
		do {
			perim += sqrt(CGAL::squared_distance(h_it->vertex()->point(), h_it->opposite()->vertex()->point()));
			++h_it;
		} while (h_it != f_it->facet_begin());
		if (perim > maxMeasure) maxMeasure = perim;
		if (perim < minMeasure) minMeasure = perim;
		map[f_it] = perim;
		moyMeasure += perim;
		++nbFacets;
		++f_it;
	}
	moyMeasure /= (double) nbFacets;
}

double computeAngle(const Point_3 & p1, const Point_3 & p2, const Point_3 & p3) {
	auto v1 = (p1 - p2);
	v1 = v1 / sqrt(v1 * v1);
	auto v2 = p3 - p2;
	v2 = v2 / sqrt(v2 * v2);
	auto n = CGAL::cross_product(v1, v2);
	return asin(n * n);
}

void calculPlusPetitAngle(Polyhedron &mesh, Facet_double_map &map) {
	Facet_iterator f_it = mesh.facets_begin();
	std::vector<Point_3> points;
	double angle1, angle2, angle3, minAngle;
	int nbFacets = 0;
	while (f_it != mesh.facets_end()) {
		minAngle = DBL_MAX;
		Halfedge_facet_circulator h_it = f_it->facet_begin();
		CGAL_assertion(CGAL::circulator_size(h_it) >= 3);
		do {
			points.push_back(h_it->vertex()->point());
			++h_it;
		} while (h_it != f_it->facet_begin());
		angle1 = computeAngle(points[0], points[1], points[2]);
		angle2 = computeAngle(points[1], points[2], points[0]);
		angle3 = computeAngle(points[2], points[0], points[1]);
		minAngle = angle1;
		if (angle2 < minAngle) minAngle = angle2;
		if (angle3 < minAngle) minAngle = angle3;
		map[f_it] = minAngle;
		if (minAngle > maxMeasure) maxMeasure = minAngle;
		if (minAngle < minMeasure) minMeasure = minAngle;
		moyMeasure += minAngle;
		++nbFacets;
		++f_it;
	}
	moyMeasure /= (double) nbFacets;
}

void calculAireOtsu(Polyhedron &mesh, Facet_double_map &map) {
	Facet_iterator f_it = mesh.facets_begin();
	Vect_double dist;
	double p, aire;
	while (f_it != mesh.facets_end()) {
		Halfedge_facet_circulator h_it = f_it->facet_begin();
		CGAL_assertion(CGAL::circulator_size(h_it) >= 3);
		do {
			dist.push_back(sqrt(CGAL::squared_distance(h_it->vertex()->point(), h_it->opposite()->vertex()->point())));
			++h_it;
		} while (h_it != f_it->facet_begin());
		p = (dist[0] + dist[1] + dist[2])/2;
		aire = sqrt(p * (p - dist[0]) * (p - dist[1]) * (p - dist[2]));
		map[f_it] = aire;
		if (aire > maxArea) maxArea = aire;
		if (aire < minArea) minArea = aire;
		++f_it;
	}
}

Vect_double randColor() {
	Vect_double color;
	for (int i = 0 ; i < 3 ; ++i) color.push_back((double) std::rand() / (double) RAND_MAX);
	return color;
}

void calculClasses(Polyhedron &mesh, Facet_double_map &measure_map, Facet_int_map &color_map) {
	Facet_iterator f_it = mesh.facets_begin();
	while (f_it != mesh.facets_end()) {
		if (measure_map[f_it] < moyMeasure) color_map[f_it] = 0;
		else color_map[f_it] = 1;
		++f_it;
	}
}

void seuilByOtsu(Polyhedron &mesh, Facet_double_map &measure_map, Facet_int_map &color_map) {
	Vect_int histo = Vect_int(64);
	Vect_int::iterator it = histo.begin();
	int section, threshold = 0;
	int nbFacets = measure_map.size();
	int sum = 0, sumB = 0, q1 = 0, q2;
	double u1, u2, current, max = 0.0;
	Facet_double_map aire_map;
	calculAireOtsu(mesh, aire_map);
	while (it != histo.end()) {
		*it = 0;
		++it;
	}
	Facet_iterator f_it = mesh.facets_begin();
	while (f_it != mesh.facets_end()) {
		section = ((measure_map[f_it] - minMeasure) / (maxMeasure - minMeasure)) * 64;
		histo[section] += ((aire_map[f_it] - minArea) / (maxArea - minArea)) * 100;
		++f_it;
	}
	for (int i = 0; i < 63; ++i) {
		sum += i * histo[i];
	}
	for (int i = 0; i < 63; ++i) {
		q1 += histo[i];
		if (q1 == 0) continue;
		q2 = nbFacets - q1;
		sumB += i * histo[i];
		u1 = (double)sumB / (double)q1;
		u2 = (double)(sum - sumB) / (double)q2;
		current = (double)q1 * (double)q2 * (u1 - u2) * (u1 - u2);
		if (current > max) {
			threshold = i;
			max = current;
		}
	}
	for (f_it = mesh.facets_begin() ; f_it != mesh.facets_end() ; ++f_it) {
		section = (measure_map[f_it] - minMeasure) / (maxMeasure- minMeasure) * 64;
		if (section < threshold) color_map[f_it] = 0;
		else color_map[f_it] = 1;
	}
}

int segmentationComposantesConnexes(Polyhedron &mesh, Facet_int_map &map, Facet_int_map &def_map) {
	int numComposante = 0;
	int segCourante;
	Facet_Queue queue;
	Facet_bool_map parcours;
	Facet_iterator f_it = mesh.facets_begin();
	Polyhedron::Facet_handle s;
	while (f_it != mesh.facets_end()) {
		parcours[f_it++] = false;
	}
	for (f_it = mesh.facets_begin() ; f_it != mesh.facets_end() ; ++f_it) {
		if (!parcours[f_it]) {
			segCourante = map[f_it];
			queue.push(f_it);
			parcours[f_it] = true;
			while (!queue.empty()) {
				s = queue.front();
				queue.pop();
				def_map[s] = numComposante;
				Halfedge_facet_circulator h_it = s->facet_begin();
				do {
					Polyhedron::Facet_handle f = h_it->opposite()->facet();
					if (!parcours[f] && map[f] == segCourante) {
						queue.push(f);
						parcours[f] = true;
					}
					++h_it;
				} while (h_it != s->facet_begin());
			}
			numComposante++;
		}
	}
	return numComposante;
}

void colorMesh(Polyhedron &mesh, Facet_int_map &map, const char *filename, int nbComposantes) {
	std::fstream output(filename);
	output << "OFF" << std::endl; //ligne d'entête
	output << mesh.size_of_vertices() << " " << mesh.size_of_facets() << " 0" << std::endl; //infos sur le mesh
	std::copy(mesh.points_begin(), mesh.points_end(), std::ostream_iterator<Point_3>(output, "\n"));
	Vect_Color vect_color;
	for (int i = 0; i < nbComposantes; ++i) {
		vect_color.push_back(randColor());
	}
	for (Facet_iterator it = mesh.facets_begin(); it != mesh.facets_end(); ++it) {
		Halfedge_facet_circulator j = it->facet_begin();
        CGAL_assertion(CGAL::circulator_size(j) >= 3);
        output << CGAL::circulator_size(j) << ' ';
        do {
            output << ' ' << std::distance(mesh.vertices_begin(), j->vertex());
        } while ( ++j != it->facet_begin());
		output << " " << std::setprecision(2) << vect_color[map[it]][0]
			   << " " << std::setprecision(2) << vect_color[map[it]][1] 
			   << " " << std::setprecision(2) << vect_color[map[it]][2] << std::endl;
	}
}

int main(int argc, char *argv[]) {
	if (argc < 4) {
		std::cerr << "Il manque un paramètre au programme. Veuillez lui donner en entrée un nom de fichier au format off, " 
		          << "une mesure de segmentation (p pour périmètre, a pour angle) et un type de seuillage (n pour normal, o pour otsu)." << std::endl;
		return 1;
	}
	
	Polyhedron mesh;
	std::ifstream input(argv[1]);
	if (!input || !(input >> mesh) || mesh.is_empty()) {
		std::cerr << "Le fichier donné n'est pas un fichier off valide." << std::endl;
		return 1;
	}
	
	char mesure = *(argv[2]);
	if (mesure != 'p' && mesure != 'a') {
		std::cerr << "La mesure de segmentation est incorrecte." << std::endl;
		return 1;
	}
	
	char type = *(argv[3]);
	if (type != 'n' && type != 'o') {
		std::cerr << "Le type de seuillage est incorrect." << std::endl;
		return 1;
	}
	
	std::srand(std::time(0));
	
	Facet_double_map measure_map;
	Facet_int_map color_map, def_color_map;
	int nbComposantes;
	
	if (mesure == 'p') {
		calculPerimetre(mesh, measure_map);
	} else if (mesure =='a') {
		calculPlusPetitAngle(mesh, measure_map);
	}
	
	if (type == 'n') {
		calculClasses(mesh, measure_map, color_map);
	} else if (type == 'o') {
		seuilByOtsu(mesh, measure_map, color_map);
	}
	
	nbComposantes = segmentationComposantesConnexes(mesh, color_map, def_color_map);
	std::cout << "Nb Composantes : " << nbComposantes << std::endl;
	
	std::string filename = "resultMesh.off";
	colorMesh(mesh, def_color_map, filename.c_str(), nbComposantes);
	
	return 0;
}
