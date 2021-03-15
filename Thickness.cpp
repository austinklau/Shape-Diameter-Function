#include <math.h>
#include <random>
#include <numeric>
#include <iostream>
#include <ctime>
#include "Thickness.h"
#include "mesh.h"
#include "KDTree.h"
#include "KDTreeSAH.h"

using namespace std;


namespace core {

	SDF::SDF()
	{
	}

	const double threshold = 100.0;

	// Algorithm to compute SDF of a mesh at every point
	void SDF::computeSDF(Mesh& mesh, double angle, int rayCount, double epsilon, vector<double>& result) {
		int misses = 0;
		int overflow = 0;
		mesh.GenerateVertexNormal();
		mesh.GenerateFaceNormal(false, true);
		mesh.GenerateEdgeAndAdjacency();

		int numFaces = mesh.GetNumFace();
		vector<Vector3> points = mesh.GetPoint();
		vector<Vector3> faceNorms = mesh.GetFaceNormal();
		vector<Vector3> vertexNorms = mesh.GetVertexNormal();


		KDTreeSAH tree;
		tree.Build(mesh); // build KDTree
		srand(time(NULL)); // init RNG

		// Iterate over all the points
		for (int i = 0; i < points.size(); i++) {
			vector<int> neighbors;
			mesh.GetVertexStar(i, neighbors); // get all neighboring FACES of vertex i

			Vector3 rayNormal = vertexNorms[i]; // normal at vertex
			rayNormal.Negate(); // point it inwards
			rayNormal.Normalize();
			Vector3 rayOrigin = points[i]; // origin for raycast

			// save subdistances and weights here for statistical analysis
			vector<double> distances;
			vector<double> weights;

			// u and v define a plane perpindicular to rayNormal, used later for generating random rays
			Vector3 u = Vector3(-rayNormal.z(), 0, rayNormal.x());
			u.Normalize();
			Vector3 v = rayNormal.Cross(u);
			v.Normalize();

			// Cast 'rayCount' rays in a cone defined by angle 'angle'
			for (int j = 0; j < rayCount; j++) {
				// Random point generation
				double phi = randomDouble(0, 2 * PI);
				double z = randomDouble(cos(angle), 1);
				double theta = acos(z);
				Vector3 randomRay = (u * cos(phi) + v * sin(phi)) * sin(theta) + rayNormal * cos(theta); // random vector in cone
				randomRay.Normalize();

				// Do ray cast
				int faceID;
				double subDistance = 0;
				bool hit = tree.Intersect(rayOrigin, randomRay, subDistance, faceID, epsilon);

				// If successful cast & not a neighbor of the origin vertex
				if (hit && std::find(neighbors.begin(), neighbors.end(), faceID) == neighbors.end()) {
					Vector3 faceNormal = faceNorms[faceID];
					faceNormal.Normalize();
					double dp = faceNormal.Dot(randomRay); // larger = better
					double angleToOrigin = acos(faceNormal.Dot(-rayNormal));

					//// Experimental multisampling code begin (this increases runtime)
					//vector<int> faceNeighbors;
					//vector<int> vertexNeighbors;
					//mesh.GetFaceStar(faceID, faceNeighbors);
					//mesh.ExpandFaceIdIntoIndex(faceNeighbors, vertexNeighbors);
					//sort(vertexNeighbors.begin(), vertexNeighbors.end());
					//vertexNeighbors.erase(unique(vertexNeighbors.begin(), vertexNeighbors.end()), vertexNeighbors.end());

					//// Get the best distance from all vertices in the neighborhood, based on max dot product

					//for (int k : vertexNeighbors) {
					//	Vector3 point = points[k];
					//	Vector3 norm = vertexNorms[k];
					//	double subdp = norm.Dot(randomRay);
					//	if (subdp > dp) {
					//		dp = subdp;
					//		subDistance = distanceBetween(rayOrigin, point);
					//	}
					//}
					//// Experimental multisampling code end

					double weight = 1 / theta;
					if (angleToOrigin >= PI / 2 && weight < threshold) {
						distances.push_back(subDistance);
						weights.push_back(weight);
					}
				}
			}


			// Statistical Analysis
			if (distances.size() == 0) { // miss case
				result.push_back(-1);
				misses++;
				continue;
			}
			if (distances.size() == 1) {
				result.push_back(distances[0]);
				continue;
			}
			double sum = accumulate(distances.begin(), distances.end(), 0.0);
			double mean = sum / distances.size();
			double sq_sum = inner_product(distances.begin(), distances.end(), distances.begin(), 0.0);
			double stdev = sqrt(sq_sum / distances.size() - mean * mean); // standard deviation
			double median = vectorMedian(distances);

			// Compute weighted avg of non-outliers
			double total = 0.0, divisor = 0.0;

			for (int j = 0; j < distances.size(); j++) {
				if (distances[j] >= median - stdev && distances[j] <= median + stdev) {
					total += distances[j] * weights[j];
					divisor += weights[j];
				}
			}

			// Testing: Compute weighted avg between Q1 and Q3
			/*sort(distances.begin(), distances.end());
			double s = distances.size();
			for (int j = s / 4; j < 3 * s / 4; j++) {
				total += distances[j] * weights[j];
				divisor += weights[j];
			}*/

			total /= divisor;
			if (divisor == 0) total = median;


			if (total > 0 && total < threshold) result.push_back(total);
			else {
				result.push_back(-2);
				overflow++;
			}

		}

		int totalmiss = misses + overflow;
		double ratio = (1.0 - (double)totalmiss / points.size()) * 100;
		// This line can be commented out if these stats are not desired
		std::cout << misses << " Misses " << overflow << " Over Threshold | Success Ratio = " << ratio << "%" << std::endl;
	}

	// returns a random number within the interval (min, max)
	double SDF::randomDouble(double min, double max) {
		return (max - min) * ((double)rand() / (double)RAND_MAX) + min;
	}

	// median of a vector, assuming it has 2+ elements
	// (using nth_element is faster than sorting)
	double SDF::vectorMedian(vector<double>& v) {
		int n = v.size() / 2;
		if (n % 2 == 0) {
			const auto median_it1 = v.begin() + n / 2 - 1;
			const auto median_it2 = v.begin() + n / 2;

			nth_element(v.begin(), median_it1, v.end());
			const auto e1 = *median_it1;

			nth_element(v.begin(), median_it2, v.end());
			const auto e2 = *median_it2;

			return (e1 + e2) / 2;
		}
		else {
			const auto median_it = v.begin() + n / 2;
			nth_element(v.begin(), median_it, v.end());
			return *median_it;
		}

	}
	double SDF::distanceBetween(const Vector3& u, const Vector3& v) {
		double intermediate = (u.x() - v.x()) * (u.x() - v.x()) +
			(u.y() - v.y()) * (u.y() - v.y()) +
			(u.z() - v.z()) * (u.z() - v.z());
		return sqrt(intermediate);
	}
}
