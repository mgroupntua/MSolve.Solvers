using System;
using System.Collections.Generic;
using System.Linq;
using MGroup.MSolve.Discretization;

//TODO: Extend this works to 3D cells
namespace MGroup.MSolve.Meshes.Unstructured.Gmsh
{
	/// <summary>
	/// Converts cell types and the order of their vertices from GMSH to MSolve.
	/// </summary>
	internal class GmshCellFactory
	{
		private static readonly IReadOnlyDictionary<int, CellType> gmshCellCodes;

		// Vertex order for cells. Index = gmsh order, value = MSolve order.
		private static readonly IReadOnlyDictionary<CellType, int[]> gmshCellConnectivity;

		static GmshCellFactory()
		{
			var codes = new Dictionary<int, CellType>();
			codes.Add(2, CellType.Tri3);
			codes.Add(3, CellType.Quad4);
			codes.Add(9, CellType.Tri6);
			codes.Add(10, CellType.Quad9);
			codes.Add(16, CellType.Quad8);
			gmshCellCodes = codes;

			var connectivity = new Dictionary<CellType, int[]>();
			connectivity.Add(CellType.Tri3, new int[] { 0, 1, 2 });
			connectivity.Add(CellType.Quad4, new int[] { 0, 1, 2, 3 });
			connectivity.Add(CellType.Tri6, new int[] { 0, 1, 2, 3, 4, 5 });
			connectivity.Add(CellType.Quad9, new int[] { 0, 1, 2, 3, 4, 5, 6, 7, 8 });
			connectivity.Add(CellType.Quad8, new int[] { 0, 1, 2, 3, 4, 5, 6, 7 });
			gmshCellConnectivity = connectivity;
		}

		private readonly IReadOnlyList<double[]> allNodes;

		public GmshCellFactory(IReadOnlyList<double[]> allNodes)
		{
			this.allNodes = allNodes;
		}

		/// <summary>
		/// Returns true and a <see cref="CellConnectivity"/> if the <paramref name="cellCode"/> corresponds to a valid 
		/// MSolve <see cref="CellType"/>. 
		/// Otherwise returns false and null.
		/// </summary>
		/// <param name="cellCode"></param>
		/// <param name="nodeIDs"> These must be 0-based</param>
		/// <param name="cell"></param>
		/// <returns></returns>
		public int[] FindElementNodes(CellType cellType, int[] gmshNodeIDs)
		{
			int[] msolveNodeIDs = Permute(gmshNodeIDs, gmshCellConnectivity[cellType]);
			List<double[]> nodeCoords = msolveNodeIDs.Select(n => allNodes[n]).ToList();
			int[] fixedNodeOrder = FixNodeOrder(nodeCoords);
			return Permute(msolveNodeIDs, fixedNodeOrder);
		}

		public bool TryGetMSolveCellType(int gmshCellCode, out CellType cellType) 
			=> gmshCellCodes.TryGetValue(gmshCellCode, out cellType);

		/// <summary>
		/// If the order is clockwise, it is reversed.
		/// </summary>
		/// <param name="cellVertices"></param>
		private int[] FixNodeOrder(IList<double[]> cellVertices)
		{
			//TODO: This only works for 2D cells. Not sure if it sufficient or required for second order elements.
			// The area of a cell with clockwise vertices is negative.
			double cellArea = 0.0; // Actually double the area will be computed, but we only care about the sign here
			for (int i = 0; i < cellVertices.Count; ++i)
			{
				double[] vertex1 = cellVertices[i];
				double[] vertex2 = cellVertices[(i + 1) % cellVertices.Count];
				cellArea += vertex1[0] * vertex2[1] - vertex2[0] * vertex1[1];
			}

			int[] order = Enumerable.Range(0, cellVertices.Count).ToArray();
			if (cellArea < 0)
			{
				Array.Reverse(order);
			}
			return order;
		}

		private static T[] Permute<T>(T[] originalAray, int[] permutation)
		{
			var result = new T[originalAray.Length];
			for (int i = 0; i < originalAray.Length; ++i)
			{
				result[permutation[i]] = originalAray[i];
			}
			return result;
		}
	}
}
