using System;
using System.Collections.Generic;
using System.Text;
using MGroup.MSolve.Discretization;

namespace MGroup.MSolve.Meshes.Unstructured
{
	public class UnstructuredMesh : IUnstructuredMesh
	{
		public List<double[]> Nodes { get; } = new List<double[]>();

		public List<(CellType cellType, int[] nodes)> Elements { get; } = new List<(CellType, int[])>();

		public IEnumerable<(int elementID, CellType cellType, int[] nodeIDs)> EnumerateElements()
		{
			for (int e = 0; e < Elements.Count; ++e)
			{
				yield return (e, Elements[e].cellType, Elements[e].nodes);
			}
		}

		public IEnumerable<(int nodeID, double[] coordinates)> EnumerateNodes()
		{
			for (int n = 0; n < Nodes.Count; ++n)
			{
				yield return (n, Nodes[n]);
			}
		}
	}
}
