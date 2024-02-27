using System.Diagnostics;
using MGroup.MSolve.Discretization;

namespace MGroup.MSolve.Meshes.Structured
{
	public class UniformCartesianMesh3D : ICartesianMesh
	{
		private const int dim = 3;
		private const int numNodesPerElement = 8;

		private readonly int axisMajor;
		private readonly int axisMedium;
		private readonly int axisMinor;
		private readonly double[] dx;
		private readonly int[][] elementNodeIdxOffsets;
		private readonly int firstNodeID;
		private readonly int firstElementID;

		private UniformCartesianMesh3D(double[] minCoordinates, double[] maxCoordinates, int[] numElements, 
			int axisMajor, int axisMedium, int axisMinor, int[] elementNodeOrderPermutation, int firstNodeID, int firstElementID)
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

			NumNodesTotal = NumNodes[0] * NumNodes[1] * NumNodes[2];
			NumElementsTotal = NumElements[0] * NumElements[1] * NumElements[2];

			this.axisMajor = axisMajor;
			this.axisMedium = axisMedium;
			this.axisMinor = axisMinor;
			this.firstNodeID = firstNodeID;
			this.firstElementID = firstElementID;

			// Hexa8 node order
			var elementNodeIdxOffsetsDefault = new int[8][];
			elementNodeIdxOffsetsDefault[0] = new int[] { 0, 0, 0 };
			elementNodeIdxOffsetsDefault[1] = new int[] { 1, 0, 0 };
			elementNodeIdxOffsetsDefault[2] = new int[] { 1, 1, 0 };
			elementNodeIdxOffsetsDefault[3] = new int[] { 0, 1, 0 };
			elementNodeIdxOffsetsDefault[4] = new int[] { 0, 0, 1 };
			elementNodeIdxOffsetsDefault[5] = new int[] { 1, 0, 1 };
			elementNodeIdxOffsetsDefault[6] = new int[] { 1, 1, 1 };
			elementNodeIdxOffsetsDefault[7] = new int[] { 0, 1, 1 };

			this.elementNodeIdxOffsets = new int[numNodesPerElement][];
			for (int n = 0; n < numNodesPerElement; ++n)
			{
				this.elementNodeIdxOffsets[elementNodeOrderPermutation[n]] = elementNodeIdxOffsetsDefault[n];
			}
		}

		public CellType CellType => CellType.Hexa8;

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
			for (int k = 0; k < NumNodes[axisMinor]; ++k)
			{
				for (int j = 0; j < NumNodes[axisMedium]; ++j)
				{
					for (int i = 0; i < NumNodes[axisMajor]; ++i)
					{
						var idx = new int[dim];
						idx[axisMinor] = k;
						idx[axisMedium] = j;
						idx[axisMajor] = i;

						yield return (GetNodeID(idx), GetNodeCoordinates(idx));
					}
				}
			}
		}

		public IEnumerable<(int elementID, int[] nodeIDs)> EnumerateElements()
		{
			for (int k = 0; k < NumElements[axisMinor]; ++k)
			{
				for (int j = 0; j < NumElements[axisMedium]; ++j)
				{
					for (int i = 0; i < NumElements[axisMajor]; ++i)
					{
						var idx = new int[dim];
						idx[axisMinor] = k;
						idx[axisMedium] = j;
						idx[axisMajor] = i;

						yield return (GetElementID(idx), GetElementConnectivity(idx));
					}
				}
			}
		}

		public int GetNodeID(int[] nodeIdx)
		{
			CheckNodeIdx(nodeIdx);

			// E.g. x-major, y-medium, z-minor: id = iX + iY * numNodesX + iZ * NumNodesX * NumNodesY
			return firstNodeID + nodeIdx[axisMajor] + nodeIdx[axisMedium] * NumNodes[axisMajor] 
				+ nodeIdx[axisMinor] * NumNodes[axisMajor] * NumNodes[axisMedium];
		}

		public int[] GetNodeIdx(int nodeID)
		{
			CheckNodeID(nodeID);

			int id = nodeID - firstNodeID;
			int numNodesPlane = NumNodes[axisMajor] * NumNodes[axisMedium];
			int mod = id % numNodesPlane;

			var idx = new int[dim];
			idx[axisMinor] = id / numNodesPlane;
			idx[axisMedium] = mod / NumNodes[axisMajor];
			idx[axisMajor] = mod % NumNodes[axisMajor];

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

			// E.g. x-major, y-medium, z-minor: id = iX + iY * NumElementsX + iZ * NumElementsX * NumElementsY
			return firstElementID + elementIdx[axisMajor] + elementIdx[axisMedium] * NumElements[axisMajor]
				+ elementIdx[axisMinor] * NumElements[axisMajor] * NumElements[axisMedium];
		}

		public int[] GetElementIdx(int elementID)
		{
			CheckElementID(elementID);

			int id = elementID - firstElementID;
			int numElementsPlane = NumElements[axisMajor] * NumElements[axisMedium];
			int mod = id % numElementsPlane;

			var idx = new int[dim];
			idx[axisMinor] = id / numElementsPlane;
			idx[axisMedium] = mod / NumElements[axisMajor];
			idx[axisMajor] = mod % NumElements[axisMajor];

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
				nodeIdx[2] = elementIdx[2] + offset[2];
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
			private int axisMinorChoice;
			private int[] elementNodeOrderPermutation;
			private int firstElementID;
			private int firstNodeID;
			private bool isAxisMajorMinorDefault;

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
				isAxisMajorMinorDefault = true;
				axisMajorChoice = int.MinValue;
				axisMinorChoice = int.MinValue;
				elementNodeOrderPermutation = new int[] { 0, 1, 2, 3, 4, 5, 6, 7 };
				firstNodeID = 0;
				firstElementID = 0;
			}

			public UniformCartesianMesh3D BuildMesh()
			{
				Validate();
				(int majorAxis, int mediumAxis, int minorAxis) = ChooseAxes();
				ValidateAxes(majorAxis, mediumAxis, minorAxis);

				return new UniformCartesianMesh3D(coordsMin.Copy(), coordsMax.Copy(), numElements.Copy(),
					majorAxis, mediumAxis, minorAxis, elementNodeOrderPermutation.Copy(), firstNodeID, firstElementID);
			}

			/// <summary>
			/// Applies Bathe's order of nodes in each Hexa8 element of the mesh:
			/// node0 = (+1, +1, +1), node1 = (-1, +1, +1), node2 = (-1, -1, +1), node3 = (+1, -1, +1), 
			/// node4 = (+1, +1, -1), node5 = (-1, +1, -1), node6 = (-1, -1, -1), node7 = (+1, -1, -1) 
			/// where the +-1 coordinates correspond to the local coordinate system of the element.
			/// </summary>
			/// <returns>This object for chaining.</returns>
			public Builder SetElementNodeOrderBathe()
			{
				this.elementNodeOrderPermutation = new int[] { 6, 7, 4, 5, 2, 3, 0, 1 };
				return this;
			}

			/// <summary>
			/// Reapplies the default order of nodes in each Hexa8 element of the mesh:
			/// node0 = (-1, -1, -1), node1 = (+1, -1, -1), node2 = (+1, +1, -1), node3 = (-1, +1, -1), 
			/// node4 = (-1, -1, +1), node5 = (+1, -1, +1), node6 = (+1, +1, +1), node7 = (-1, +1, +1)
			/// where the +-1 coordinates correspond to the local coordinate system of the element.
			/// </summary>
			/// <returns>This object for chaining.</returns>
			public Builder SetElementNodeOrderDefault()
			{
				this.elementNodeOrderPermutation = new int[] { 0, 1, 2, 3, 4, 5, 6, 7 };
				return this;
			}

			/// <summary>
			/// Configures the order of nodes in each Quad4 element of the mesh. By default this order is: 
			/// node0 = (-1, -1, -1), node1 = (+1, -1, -1), node2 = (+1, +1, -1), node3 = (-1, +1, -1), 
			/// node4 = (-1, -1, +1), node5 = (+1, -1, +1), node6 = (+1, +1, +1), node7 = (-1, +1, +1)
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
			/// The node ids will be ordered such that they are contiguous along dimension <paramref name="majorAxis"/>, while they  
			/// will have the maximum id difference along dimension <paramref name="minorAxis"/>. 
			/// Calling this method overrides the default node order: nodes are contiguous in the dimension with mininum
			/// number of nodes and have the maximum id difference in the dimension with the maximum number of nodes.
			/// </summary>
			/// <param name="majorAxis">
			/// The axis along which node ids will be contiguous. 0 for x, 1 for y or 2 for z. 
			/// Must be different from <paramref name="minorAxis"/>.
			/// </param>
			/// <param name="minorAxis">
			/// The axis along which node ids will have the maximum distance (compared to other axes). 
			/// 0 for x, 1 for y or 2 for z. Must be different from <paramref name="majorAxis"/>.
			/// </param>
			/// <returns>This object for chaining.</returns>
			public Builder SetMajorMinorAxis(int majorAxis, int minorAxis)
			{
				this.isAxisMajorMinorDefault = false;
				this.axisMajorChoice = majorAxis;
				this.axisMinorChoice = minorAxis;
				return this;
			}

			/// <summary>
			/// Reapplies the default order of node ids: nodes are contiguous in the dimension with mininum number of nodes
			/// and have the maximum id difference in the dimension with the maximum number of nodes.
			/// </summary>
			/// <returns>This object for chaining.</returns>
			public Builder SetMajorMinorAxisDefault()
			{
				this.isAxisMajorMinorDefault = true;
				this.axisMajorChoice = int.MinValue;
				this.axisMinorChoice = int.MinValue;
				return this;
			}

			private (int majorAxis, int mediumAxis, int minorAxis) ChooseAxes()
			{
				int majorAxis, mediumAxis, minorAxis;
				if (isAxisMajorMinorDefault)
				{
					// Decide based on the number of nodes/elements per axis 
					// Sort axes based on their number of elements
					var entries = new List<(int count, int axis)>();
					entries.Add((numElements[0], 0));
					entries.Add((numElements[1], 1));
					entries.Add((numElements[2], 2));
					int[] sortedAxes = SortAxes(entries);

					majorAxis = sortedAxes[0];
					mediumAxis = sortedAxes[1];
					minorAxis = sortedAxes[2];
				}
				else
				{
					// Respect client decision
					majorAxis = axisMajorChoice; 
					minorAxis = axisMinorChoice;
					if ((majorAxis == 0) && (minorAxis == 1)) mediumAxis = 2;
					else if ((majorAxis == 0) && (minorAxis == 2)) mediumAxis = 1;
					else if ((majorAxis == 1) && (minorAxis == 0)) mediumAxis = 2;
					else if ((majorAxis == 1) && (minorAxis == 2)) mediumAxis = 0;
					else if ((majorAxis == 2) && (minorAxis == 0)) mediumAxis = 1;
					else if ((majorAxis == 2) && (minorAxis == 1)) mediumAxis = 0;
					else throw new ArgumentException("Major and minors axes must be 0, 1 or 2 and different from each other");
				}
				return (majorAxis, mediumAxis, minorAxis);
			}

			private static int[] SortAxes(List<(int count, int axis)> entries)
			{
				Debug.Assert(entries.Count == 3);
				var sortedAxes = new int[3];
				int idx = 0;
				while (idx < 3)
				{
					int min = int.MaxValue;
					int axisOfMin = -1;

					foreach ((int count, int axis) in entries)
					{
						if (count < min)
						{
							min = count;
							axisOfMin = axis;
						}
						else if (count == min) // prefer axis x over y, z and axis y over z
						{
							if (axis < axisOfMin)
							{
								min = count;
								axisOfMin = axis;
							}
						}
					}

					sortedAxes[idx++] = axisOfMin;
					bool removedEntry = entries.Remove((min, axisOfMin));
					Debug.Assert(removedEntry);
				}
				return sortedAxes;
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

				// Major, medium, minor axes will be checked in a different method
				// First node & element IDs: no illegal values yet.
			}

			private static void ValidateAxes(int majorAxis, int mediumAxis, int minorAxis)
			{
				bool isCorrect = false;
				isCorrect |= (majorAxis == 0) && (minorAxis == 1) && (mediumAxis == 2);
				isCorrect |= (majorAxis == 0) && (minorAxis == 2) && (mediumAxis == 1);
				isCorrect |= (majorAxis == 1) && (minorAxis == 0) && (mediumAxis == 2);
				isCorrect |= (majorAxis == 1) && (minorAxis == 2) && (mediumAxis == 0);
				isCorrect |= (majorAxis == 2) && (minorAxis == 0) && (mediumAxis == 1);
				isCorrect |= (majorAxis == 2) && (minorAxis == 1) && (mediumAxis == 0);

				if (!isCorrect)
				{
					throw new ArgumentException("Major and minors axes must be 0, 1 or 2 and different from each other");
				}
			}
		}
	}
}
