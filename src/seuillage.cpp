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
typedef Kernel::Point_3 Point_3;
typedef Polyhedron::Facet_iterator Facet_iterator;
typedef Polyhedron::Vertex_iterator Vertex_iterator;
typedef Polyhedron::Halfedge_iterator Halfedge_iterator;
typedef Polyhedron::Halfedge_around_facet_circulator Halfedge_facet_circulator;
typedef std::map<Polyhedron::Facet_handle, double> Perim_map;
typedef std::map<Polyhedron::Facet_handle, int> Color_map;
typedef std::map<Polyhedron::Facet_handle, bool> Parcours_map;
typedef std::vector<double> Color;
typedef std::vector<Color> Vect_Color;
typedef std::queue<Polyhedron::Facet_handle> Facet_Queue;

double maxPerim = 0.0;
double minPerim = DBL_MAX;
double moyPerim = 0.0;

void calculPerimetre(Polyhedron &mesh, Perim_map &map) {
	Facet_iterator f_it = mesh.facets_begin();
	int nbFacets = 0;
	while (f_it != mesh.facets_end()) {
		Halfedge_facet_circulator h_it = f_it->facet_begin();
		CGAL_assertion(CGAL::circulator_size(h_it) >= 3);
		double perim = 0;
		do {
			perim += sqrt(CGAL::squared_distance(h_it->vertex()->point(), h_it->opposite()->vertex()->point()));
			++h_it;
		} while (h_it != f_it->facet_begin());
		if (perim > maxPerim) maxPerim = perim;
		if (perim < minPerim) minPerim = perim;
		map[f_it] = perim;
		moyPerim += perim;
		++nbFacets;
		++f_it;
	}
	moyPerim /= (double) nbFacets;
}

Color randColor() {
	Color color;
	for (int i = 0 ; i < 3 ; ++i) color.push_back((double) std::rand() / (double) RAND_MAX);
	return color;
}

void calculClasses(Polyhedron &mesh, Perim_map &perim_map, Color_map &color_map) {
	Facet_iterator f_it = mesh.facets_begin();
	while (f_it != mesh.facets_end()) {
		if (perim_map[f_it] < moyPerim) color_map[f_it] = 0;
		else color_map[f_it] = 1;
		++f_it;
	}
}

int seuilByOtsu(Polyhedron &mesh, Perim_map &perim_map, Color_map &color_map) {
	return 0;
}

int segmentationComposantesConnexes(Polyhedron &mesh, Color_map &map, Color_map &def_map) {
	int numComposante = 0;
	int segCourante;
	Facet_Queue queue;
	Parcours_map parcours;
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

void colorMesh(Polyhedron &mesh, Color_map &map, const char *filename, int nbComposantes) {
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
	if (argc < 2) {
		std::cerr << "Il manque un paramètre au programme. Veuillez lui donner en entrée un nom de fichier au format off." << std::endl;
		return 1;
	}
	
	Polyhedron mesh;
	std::ifstream input(argv[1]);
	if (!input || !(input >> mesh) || mesh.is_empty()) {
		std::cerr << "Le fichier donné n'est pas un fichier off valide." << std::endl;
		return 1;
	}
	
	std::srand(std::time(0));
	
	Perim_map perim_map;
	Color_map color_map, def_color_map;
	int nbComposantes;
	calculPerimetre(mesh, perim_map);
	calculClasses(mesh, perim_map, color_map);
	nbComposantes = segmentationComposantesConnexes(mesh, color_map, def_color_map);
	std::cout << nbComposantes << std::endl;
	
	std::string filename = "resultMesh.off";
	colorMesh(mesh, def_color_map, filename.c_str(), nbComposantes);
	
	return 0;
}
