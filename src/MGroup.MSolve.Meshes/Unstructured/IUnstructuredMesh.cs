using System;
using System.Collections.Generic;
using System.Text;
using MGroup.MSolve.Discretization;

namespace MGroup.MSolve.Meshes.Unstructured
{
	public interface IUnstructuredMesh
	{
		/// <summary>
		/// Enumerates the global IDs and coordinates of the nodes. The nodes are returned in increasing order of their id.
		/// </summary>
		IEnumerable<(int nodeID, double[] coordinates)> EnumerateNodes();

		/// <summary>
		/// For each element returns its global ID, cell type and the global IDs of its nodes. The elements are returned in 
		/// increasing order of their id.
		/// </summary>
		IEnumerable<(int elementID, CellType cellType, int[] nodeIDs)> EnumerateElements();
	}
}
