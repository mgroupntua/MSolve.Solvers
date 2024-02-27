using System.Diagnostics;

using MGroup.MSolve.Discretization;

//TODO: Allow client to control the major/minor/etc axes for numbering elements independently of vertices.
//TODO: Builder should be in a different file (and not nested). Checking the params must be done in the public constructor.
//      Invalid builder properties that will not be checked by the constructor, must be checked by the builder.
namespace MGroup.MSolve.Meshes.Structured
{
	public class UniformCartesianMesh2D : ICartesianMesh
	{
		private const int dim = 2;
		private const int numNodesPerElement = 4;

		private readonly int axisMajor;
		private readonly int axisMinor;
		private readonly double[] dx;
		private readonly int[][] elementNodeIdxOffsets;
		private readonly int firstElementID;
		private readonly int firstNodeID;

		private UniformCartesianMesh2D(double[] minCoordinates, double[] maxCoordinates, int[] numElements, int axisMajor, int axisMinor,
			int[] elementNodeOrderPermutation, int firstNodeID, int firstElementID)
		{
			this.MinCoordinates = minCoordinates;
			this.MaxCoordinates = maxCoordinates;
			this.NumElements = numElements;

			NumNodes = new int[dim];
			for (int d = 0; d < dim; d++)
			{
				NumNodes[d] = numElements[d] + 1;
			}

			dx = new double[dim];
			for (int d = 0; d < dim; d++)
			{
				dx[d] = (maxCoordinates[d] - minCoordinates[d]) / numElements[d];
			}

			NumNodesTotal = NumNodes[0] * NumNodes[1];
			NumElementsTotal = NumElements[0] * NumElements[1];

			this.axisMajor = axisMajor;
			this.axisMinor = axisMinor;
			this.firstNodeID = firstNodeID;
			this.firstElementID = firstElementID;

			// Quad4 node order
			var elementNodeIdxOffsetsDefault = new int[numNodesPerElement][];
			elementNodeIdxOffsetsDefault[0] = new int[] { 0, 0 };
			elementNodeIdxOffsetsDefault[1] = new int[] { 1, 0 };
			elementNodeIdxOffsetsDefault[2] = new int[] { 1, 1 };
			elementNodeIdxOffsetsDefault[3] = new int[] { 0, 1 };

			this.elementNodeIdxOffsets = new int[numNodesPerElement][];
			for (int n = 0; n < numNodesPerElement; ++n)
			{
				this.elementNodeIdxOffsets[elementNodeOrderPermutation[n]] = elementNodeIdxOffsetsDefault[n]; 
			}
		}

		public CellType CellType => CellType.Quad4;

		public int Dimension => dim;

		public double[] DistancesBetweenPoints => dx;

		public double[] MinCoordinates { get; }

		public double[] MaxCoordinates { get; }

		public int[] NumElements { get; }

		public int NumElementsTotal { get; }

		public int[] NumNodes { get; }

		public int NumNodesPerElement => numNodesPerElement;

		public int NumNodesTotal { get; }

		public IEnumerable<(int nodeID, double[] coordinates)> EnumerateNodes()
		{
			for (int j = 0; j < NumNodes[axisMinor]; ++j)
			{
				for (int i = 0; i < NumNodes[axisMajor]; ++i)
				{
					var idx = new int[dim];
					idx[axisMinor] = j;
					idx[axisMajor] = i;
					yield return (GetNodeID(idx), GetNodeCoordinates(idx));
				}
			}
		}

		public IEnumerable<(int elementID, int[] nodeIDs)> EnumerateElements()
		{
			for (int j = 0; j < NumElements[axisMinor]; ++j)
			{
				for (int i = 0; i < NumElements[axisMajor]; ++i)
				{
					var idx = new int[dim];
					idx[axisMinor] = j;
					idx[axisMajor] = i;
					yield return (GetElementID(idx), GetElementConnectivity(idx));
				}
			}
		}

		public int GetNodeID(int[] nodeIdx)
		{
			CheckNodeIdx(nodeIdx);

			// E.g. x-major (nodes contiguous along x): id = iX + iY * numNodesX
			return firstNodeID + nodeIdx[axisMajor] + nodeIdx[axisMinor] * NumNodes[axisMajor];
		}

		public int[] GetNodeIdx(int nodeID)
		{
			CheckNodeID(nodeID);

			// E.g. x-major (nodes contiguous along x): iX = id % numNodesX; y = id / numNodesX;
			int id = nodeID - firstNodeID;
			var idx = new int[dim];
			idx[axisMajor] = id % NumNodes[axisMajor];
			idx[axisMinor] = id / NumNodes[axisMajor];
			return idx;
		}

		public double[] GetNodeCoordinates(int[] nodeIdx)
		{
			CheckNodeIdx(nodeIdx);
			var coords = new double[dim];
			for (int d = 0; d < dim; d++)
			{
				coords[d] = MinCoordinates[d] + nodeIdx[d] * dx[d];
			}
			return coords;
		}

		public int GetElementID(int[] elementIdx)
		{
			CheckElementIdx(elementIdx);

			// E.g. x-major (elements contiguous along x): id = iX + iY * NumElementsX
			return firstElementID + elementIdx[axisMajor] + elementIdx[axisMinor] * NumElements[axisMajor];
		}

		public int[] GetElementIdx(int elementID)
		{
			CheckElementID(elementID);

			// E.g. x-major (elements contiguous along x): iX = id % numNodesX; y = id / NumElementsX;
			int id = elementID - firstElementID;
			var idx = new int[dim];
			idx[axisMajor] = id % NumElements[axisMajor];
			idx[axisMinor] = id / NumElements[axisMajor];
			return idx;
		}

		public int[] GetElementConnectivity(int[] elementIdx)
		{
			CheckElementIdx(elementIdx);
			var nodeIDs = new int[numNodesPerElement];
			var nodeIdx = new int[dim]; // Avoid allocating an array per node
			for (int n = 0; n < numNodesPerElement; ++n)
			{
				int[] offset = elementNodeIdxOffsets[n];
				nodeIdx[0] = elementIdx[0] + offset[0];
				nodeIdx[1] = elementIdx[1] + offset[1];
				nodeIDs[n] = GetNodeID(nodeIdx);
			}

			return nodeIDs;
		}

		public int[] GetElementConnectivity(int elementID) => GetElementConnectivity(GetElementIdx(elementID));

		[Conditional("DEBUG")]
		private void CheckElementIdx(int[] elementIdx)
		{
			if (elementIdx.Length != dim)
			{
				throw new ArgumentException($"Element index must be an array with Length = {dim}");
			}
			for (int d = 0; d < dim; ++d)
			{
				if ((elementIdx[d] < 0) || (elementIdx[d] >= NumElements[d]))
				{
					throw new ArgumentException($"Element index along dimension {d} must belong in [0, {NumElements[d]})");
				}
			}
		}

		[Conditional("DEBUG")]
		private void CheckElementID(int elementID)
		{
			if ((elementID < firstElementID) || (elementID >= firstElementID + NumElementsTotal))
			{
				throw new ArgumentException(
					$"Element ID must belong in [{firstElementID}, {firstElementID + NumElementsTotal})");
			}
		}

		[Conditional("DEBUG")]
		private void CheckNodeIdx(int[] nodeIdx)
		{
			if (nodeIdx.Length != dim)
			{
				throw new ArgumentException($"Node index must be an array with Length = {dim}");
			}
			for (int d = 0; d < dim; ++d)
			{
				if ((nodeIdx[d] < 0) || (nodeIdx[d] >= NumNodes[d]))
				{
					throw new ArgumentException($"Node index along dimension {d} must belong in [0, {NumNodes[d]})");
				}
			}
		}

		[Conditional("DEBUG")]
		private void CheckNodeID(int nodeID)
		{
			if ((nodeID < firstNodeID) || (nodeID >= firstNodeID + NumNodesTotal))
			{
				throw new ArgumentException(
					$"Node ID must belong in [{firstNodeID}, {firstNodeID + NumNodesTotal})");
			}
		}

		public class Builder
		{
			private readonly double[] coordsMin;
			private readonly double[] coordsMax;
			private readonly int[] numElements;
			private int axisMajorChoice;
			private int[] elementNodeOrderPermutation;
			private int firstElementID;
			private int firstNodeID;

			/// <summary>
			/// 
			/// </summary>
			/// <param name="minCoordinates"></param>
			/// <param name="maxCoordinates"></param>
			/// <param name="numElements">Array with 3 positive integers.</param>
			public Builder(double[] minCoordinates, double[] maxCoordinates, int[] numElements) 
			{
				this.coordsMin = minCoordinates.Copy();
				this.coordsMax = maxCoordinates.Copy();
				this.numElements = numElements.Copy();

				// Defaults
				axisMajorChoice = int.MinValue;
				elementNodeOrderPermutation = new int[] { 0, 1, 2, 3 };
				firstNodeID = 0;
				firstElementID = 0;
			}

			public UniformCartesianMesh2D BuildMesh()
			{
				Validate();
				(int majorAxis, int minorAxis) = ChooseAxes();
				return new UniformCartesianMesh2D(coordsMin.Copy(), coordsMax.Copy(), numElements.Copy(),
					majorAxis, minorAxis, elementNodeOrderPermutation.Copy(), firstNodeID, firstElementID);
			}

			/// <summary>
			/// Applies a clockwise order of nodes in each Quad4 element of the mesh:
			/// node0 = (-1, +1), node1 = (+1, +1), node2 = (+1, -1), node3 = (-1, -1), 
			/// where the +-1 coordinates correspond to the local coordinate system of the element.
			/// </summary>
			/// <returns>This object for chaining.</returns>
			public Builder SetElementNodeOrderClockwise()
			{
				this.elementNodeOrderPermutation = new int[] { 3, 2, 1, 0 };
				return this;
			}

			/// <summary>
			/// Reapplies the default counter-clockwise order of nodes in each Quad4 element of the mesh:
			/// node0 = (-1, -1), node1 = (+1, -1), node2 = (+1, +1), node3 = (-1, +1), 
			/// where the +-1 coordinates correspond to the local coordinate system of the element.
			/// </summary>
			/// <returns>This object for chaining.</returns>
			public Builder SetElementNodeOrderCounterClockwise()
			{
				this.elementNodeOrderPermutation = new int[] { 0, 1, 2, 3 };
				return this;
			}

			
			/// <summary>
			/// Configures the order of nodes in each Quad4 element of the mesh. By default this order is: 
			/// node0 = (-1, -1), node1 = (+1, -1), node2 = (+1, +1), node3 = (-1, +1), 
			/// where the +-1 coordinates correspond to the  local coordinate system of the element. 
			/// To change this order, <paramref name="permutation"/> must be a permutation array such that: 
			/// finalOrder[<paramref name="permutation"/>[i]] = defaultOrder[i], i = 0,1,...3 .
			/// </summary>
			/// <param name="permutation">
			/// The permutation of element nodes: finalOrder[<paramref name="permutation"/>[i]] = defaultOrder[i], i = 0,1,...3.
			/// </param>
			/// <returns>This object for chaining.</returns>
			public Builder SetElementNodeOrderPermutation(int[] permutation)
			{
				this.elementNodeOrderPermutation = permutation;
				return this;
			}

			public Builder SetFirstElementID(int elementID)
			{
				this.firstElementID = elementID;
				return this;
			}

			public Builder SetFirstNodeID(int nodeID)
			{
				this.firstNodeID = nodeID;
				return this;
			}

			/// <summary>
			/// The node IDs will be ordered such that they are contiguous along dimension <paramref name="axis"/>. Calling this
			/// method overrides the default node order: nodes are contiguous in the dimension with mininum number of nodes.
			/// </summary>
			/// <param name="axis">The axis along which node ids will be contiguous. 0 for x or 1 for y.</param>
			/// <returns>This object for chaining.</returns>
			public Builder SetMajorAxis(int axis)
			{
				this.axisMajorChoice = axis;
				return this;
			}

			/// <summary>
			/// Reapplies the default order of node ids: nodes are contiguous in the dimension with mininum number of nodes.
			/// </summary>
			/// <returns>This object for chaining.</returns>
			public Builder SetMajorAxisDefault()
			{
				this.axisMajorChoice = int.MinValue;
				return this;
			}

			private (int majorAxis, int minorAxis) ChooseAxes()
			{
				int majorAxisFinal, minorAxisFinal;
				if (axisMajorChoice == int.MinValue)
				{
					// Decide based on the number of nodes/elements per axis 
					if (numElements[0] <= numElements[1])
					{
						majorAxisFinal = 0;
						minorAxisFinal = 1;
					}
					else
					{
						majorAxisFinal = 1;
						minorAxisFinal = 0;
					}
				}
				else
				{
					// Respect client decision
					majorAxisFinal = axisMajorChoice;
					minorAxisFinal = majorAxisFinal == 0 ? 1 : 0;
				}
				return (majorAxisFinal, minorAxisFinal);
			}

			//TODO: This must be called lazily but only once when building several products with the same builder without
			//      mutating the builder.
			private void Validate()
			{
				// Coordinates
				if (coordsMin.Length != dim) throw new ArgumentException($"Length of min coordinates must be {dim}.");
				if (coordsMax.Length != dim) throw new ArgumentException($"Length of max coordinates must be {dim}.");
				for (int d = 0; d < dim; ++d)
				{
					if (coordsMin[d] >= coordsMax[d])
					{
						throw new ArgumentException(
							$"Along axis {d}, min coordinates must be strictly less than max coordinates.");
					}
				}

				// Elements per axis
				if (numElements.Length != dim) throw new ArgumentException($"Length of number of elements must be {dim}.");
				for (int d = 0; d < dim; ++d)
				{
					if (numElements[d] < 1)
					{
						throw new ArgumentException($"Along axis {d}, there must be at least 1 element.");
					}
				}

				// Major axis
				if (axisMajorChoice != int.MinValue)
				{
					if ((axisMajorChoice != 0) && (axisMajorChoice != 1)) throw new ArgumentException("Major axis must be 0 (x) or 1 (y)");
				}

				// Node order
				if (elementNodeOrderPermutation.Length != numNodesPerElement)
				{
					throw new ArgumentException("The provided permutation array for the order of nodes in each element" +
						$" must have {numNodesPerElement} entries");
				}
				if (!elementNodeOrderPermutation.AreContiguousUniqueIndices())
				{
					throw new ArgumentException("Invalid permutation array for the order of nodes in each element." +
						$" The entries must be unique and belong to [0, {numNodesPerElement})");
				}

				// First node & element IDs: no illegal values yet.
			}
		}
	}
}
