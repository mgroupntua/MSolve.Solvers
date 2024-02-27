using System;
using System.Collections.Generic;
using System.Text;

//TODO: Perhaps a similar interface for triangular structured meshes. The generalization of triangle is called simplex so perhaps 
//      ISimplicialStructuredMesh. There is actually a lot of math behind these: https://en.wikipedia.org/wiki/Simplex.
namespace MGroup.MSolve.Meshes.Structured
{
	/// <summary>
	/// Special case of structured mesh, where elements are quadrilaterals in 2D or hexahedrals in 3D.
	/// </summary>
	public interface ICartesianMesh : IStructuredMesh
	{
		double[] DistancesBetweenPoints { get; }

		int[] NumElements { get; }

		public (int elementID, double[] naturalCoords) FindElementContaining(double[] globalCoords)
		{
			int dim = Dimension;
			var elementIdx = new int[dim];
			var naturalCoords = new double[dim];
			for (int d = 0; d < dim; ++d)
			{
				// For a given axis:
				// Let x be the global cartesian coordinate of a point along this axis.
				// Let xMin be the lower bound of x.
				// Let e be the index of the element that contains x along this axis.
				// Let xi be the isoparametric coordinate of element e, that corresponds to x.
				// Let dx be the length of each element along this axis.
				// Let xStart = xMin + e*dx and xEnd = xMin + (e+1)*dx be the global cartesian coordinates of the start 
				// and end point of element e.
				// Then x = xMin + xStart * 1/2 * (1-xi) + xEnd * 1/2 * (1+xi) = ... = xMin + e*dx + dx * (1+xi)/2 
				// <=> x = xMin + xStart + dx * (1+xi)/2 <=> xi = 2*(x-xMin)/dx -2*e -1 
				// where dx * (1+xi)/2 belongs in [0, dx]
				double x = globalCoords[d];
				double xMin = MinCoordinates[d];
				double dx = DistancesBetweenPoints[d];
				elementIdx[d] = (int)Math.Floor((x - xMin) / dx);
				
				if (elementIdx[d] > NumElements[d])
				{
					throw new Exception($"Point ({globalCoords[0]}, {globalCoords[1]}, {globalCoords[2]}) is outside the mesh");
				}
				else if (elementIdx[d] == NumElements[d])
				{
					elementIdx[d] -= 1;
					naturalCoords[d] = 1.0;
				}
				else
				{
					naturalCoords[d] = 2 * (x - xMin) / dx - 2 * elementIdx[d] - 1;
				}
			}

			return (GetElementID(elementIdx), naturalCoords);
		}
	}
}
