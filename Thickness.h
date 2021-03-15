#pragma once


#include <vector>
#include "common.h"

namespace core
{
	class Mesh;
	class SDF
	{
	public:
		SDF();
		void computeSDF(Mesh& mesh, double angle, int rayCount, double epsilon, std::vector<double>& result);
	private:
		double randomDouble(double min, double max);
		double vectorMedian(std::vector<double> &v);
		double distanceBetween(const Vector3 &u, const Vector3 &v);
	
	};
	
}