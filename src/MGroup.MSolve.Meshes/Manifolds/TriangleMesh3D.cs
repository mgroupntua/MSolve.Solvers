using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Text;

namespace MGroup.MSolve.Meshes.Manifolds
{
	public class TriangleMesh3D
	{
		public List<int[]> Cells { get; } = new List<int[]>();

		public List<double[]> Vertices { get; } = new List<double[]>();

		/// <summary>
		/// The triangles will have the same normal as the rectangle, which is n = ABxAC / |ABxAC|.
		/// </summary>
		public static TriangleMesh3D CreateForRectangle(double[] pointA, double[] pointB, double[] pointC, 
			int numVerticesAlongAB, int numVerticesAlongAC)
		{
			Debug.Assert(numVerticesAlongAB >= 2);
			Debug.Assert(numVerticesAlongAC >= 2);

			// Vertices
			var dx = new double[3]; // vector along AB 
			var dy = new double[3]; // vector along AC 
			for (int d = 0; d < 3; ++d)
			{
				dx[d] = (pointB[d] - pointA[d]) / (numVerticesAlongAB - 1);
				dy[d] = (pointC[d] - pointA[d]) / (numVerticesAlongAC - 1);
			}

			var mesh = new TriangleMesh3D();
			for (int j = 0; j < numVerticesAlongAC; ++j)
			{
				for (int i = 0; i < numVerticesAlongAB; ++i)
				{
					var vertex = new double[3];
					for (int d = 0; d < 3; ++d)
					{
						vertex[d] = pointA[d] + i * dx[d] + j * dy[d];
					}
					mesh.Vertices.Add(vertex);
				}
			}

			// Cells
			for (int j = 0; j < numVerticesAlongAC - 1; ++j)
			{
				for (int i = 0; i < numVerticesAlongAB - 1; ++i)
				{
					int firstVertex = j * numVerticesAlongAB + i;

					mesh.Cells.Add(new int[]
					{
						firstVertex,
						firstVertex + 1,
						firstVertex + numVerticesAlongAB
					});

					mesh.Cells.Add(new int[]
					{
						firstVertex + 1,
						firstVertex + numVerticesAlongAB + 1,
						firstVertex + numVerticesAlongAB
					});
				}
			}

			return mesh;
		}
	}
}
