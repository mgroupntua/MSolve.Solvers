using System.Diagnostics;
using MGroup.MSolve.Discretization;

namespace MGroup.MSolve.Meshes.Structured
{
    /// <summary>
    /// Structured mesh with uniform distances between points on a cartesian grid. The elements are Tet4 (simplices) generated
    /// by dividing each Hexa8 cell of an equivalent cartesian mesh into 4 subtetrahedra. The element indices are thus int[4], 
    /// where the first 3 entries are the index of the enclosing Hexa8 of an equivalent cartesian mesh and the third entry is 
    /// 0, 1, 2 or 3 to denote one of the 4 subtetrahedra of said Hexa8.
    /// </summary>
    public class UniformSimplicialMesh3D : IStructuredMesh
    {
        private const int dim = 3;
        private const int numNodesPerElement = 4;
        private const int numSimplicesPerCartesianCell = 6;

        private readonly int axisMajor;
        private readonly int axisMedium;
        private readonly int axisMinor;
        private readonly double[] dx;
        private readonly int firstNodeID;
        private readonly int firstElementID;
        private readonly int[] numCartesianCells;

        /// <summary>
        /// Given the index (i,j,k) of the first node of a Hexa with coordinates=(-1,-1,-1), this list contains, for each 
        /// subtetrahedron of the Hexa8, the index offsets relative to (i,j,k) of the nodes of the subtetrahedron. The offset of 
        /// node (i,j,k) is (0,0,0), the offset of (i+1, j, k+1) is (1,0,1) and so forth.
        /// </summary>
        private readonly List<int[][]> elementNodeIdxOffsets;

        private UniformSimplicialMesh3D(double[] minCoordinates, double[] maxCoordinates, int[] numNodes,
            int axisMajor, int axisMedium, int axisMinor, int firstNodeID, int firstElementID)
        {
            this.MinCoordinates = minCoordinates;
            this.MaxCoordinates = maxCoordinates;
            this.NumNodes = numNodes;
            NumNodesTotal = NumNodes[0] * NumNodes[1] * NumNodes[2];

            this.numCartesianCells = new int[dim];
            for (int d = 0; d < dim; d++)
            {
                this.numCartesianCells[d] = numNodes[d] - 1;
            }
            NumElementsTotal = numSimplicesPerCartesianCell * numCartesianCells[0] * numCartesianCells[1] * numCartesianCells[2];

            dx = new double[dim];
            for (int d = 0; d < dim; d++)
            {
                dx[d] = (maxCoordinates[d] - minCoordinates[d]) / numCartesianCells[d];
            }

            this.axisMajor = axisMajor;
            this.axisMedium = axisMedium;
            this.axisMinor = axisMinor;
            this.firstNodeID = firstNodeID;
            this.firstElementID = firstElementID;

            elementNodeIdxOffsets = new List<int[][]>(numSimplicesPerCartesianCell);

            #region Tetrahedra of prism through nodes 0,1,2,3,4,5
            // The offsets of the nodes of the first tetrahedron are
            var tet0NodeIdxOffsets = new int[numNodesPerElement][];
            tet0NodeIdxOffsets[0] = new int[] { 0, 0, 0 }; // node 0
            tet0NodeIdxOffsets[1] = new int[] { 1, 0, 0 }; // node 1
            tet0NodeIdxOffsets[2] = new int[] { 1, 1, 0 }; // node 2
            tet0NodeIdxOffsets[3] = new int[] { 0, 0, 1 }; // node 4
            elementNodeIdxOffsets.Add(tet0NodeIdxOffsets);

            // The offsets of the nodes of the second tetrahedron are
            var tet1NodeIdxOffsets = new int[numNodesPerElement][];
            tet1NodeIdxOffsets[0] = new int[] { 0, 0, 0 }; // node 0
            tet1NodeIdxOffsets[1] = new int[] { 1, 1, 0 }; // node 2
            tet1NodeIdxOffsets[2] = new int[] { 0, 1, 0 }; // node 3
            tet1NodeIdxOffsets[3] = new int[] { 0, 0, 1 }; // node 4
            elementNodeIdxOffsets.Add(tet1NodeIdxOffsets);

            // The offsets of the nodes of the third tetrahedron are
            var tet2NodeIdxOffsets = new int[numNodesPerElement][];
            tet2NodeIdxOffsets[0] = new int[] { 1, 0, 0 }; // node 1
            tet2NodeIdxOffsets[1] = new int[] { 1, 1, 0 }; // node 2
            tet2NodeIdxOffsets[2] = new int[] { 0, 0, 1 }; // node 4
            tet2NodeIdxOffsets[3] = new int[] { 1, 0, 1 }; // node 5
            elementNodeIdxOffsets.Add(tet2NodeIdxOffsets);
            #endregion

            #region Tetrahedra of prism through nodes 4,5,6,7,2,3
            // The offsets of the nodes of the fourth tetrahedron are
            var tet3NodeIdxOffsets = new int[numNodesPerElement][];
            tet3NodeIdxOffsets[0] = new int[] { 1, 1, 1 }; // node 6
            tet3NodeIdxOffsets[1] = new int[] { 1, 0, 1 }; // node 5
            tet3NodeIdxOffsets[2] = new int[] { 0, 0, 1 }; // node 4
            tet3NodeIdxOffsets[3] = new int[] { 1, 1, 0 }; // node 2
            elementNodeIdxOffsets.Add(tet3NodeIdxOffsets);

            // The offsets of the nodes of the fifth tetrahedron are
            var tet4NodeIdxOffsets = new int[numNodesPerElement][];
            tet4NodeIdxOffsets[0] = new int[] { 1, 1, 1 }; // node 6
            tet4NodeIdxOffsets[1] = new int[] { 0, 0, 1 }; // node 4
            tet4NodeIdxOffsets[2] = new int[] { 0, 1, 1 }; // node 7
            tet4NodeIdxOffsets[3] = new int[] { 1, 1, 0 }; // node 2
            elementNodeIdxOffsets.Add(tet4NodeIdxOffsets);

            // The offsets of the nodes of the sixth tetrahedron are
            var tet5NodeIdxOffsets = new int[numNodesPerElement][];
            tet5NodeIdxOffsets[0] = new int[] { 0, 1, 0 }; // node 3
            tet5NodeIdxOffsets[1] = new int[] { 0, 0, 1 }; // node 4
            tet5NodeIdxOffsets[2] = new int[] { 1, 1, 0 }; // node 2
            tet5NodeIdxOffsets[3] = new int[] { 0, 1, 1 }; // node 7
            elementNodeIdxOffsets.Add(tet5NodeIdxOffsets);
            #endregion
        }

        public CellType CellType => CellType.Tet4;

        public int Dimension => dim;

        public double[] MinCoordinates { get; }

        public double[] MaxCoordinates { get; }

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
            for (int elementID = firstElementID; elementID < firstElementID + NumElementsTotal; ++elementID)
            {
                yield return (elementID, GetElementConnectivity(elementID));
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
            int cartesianCellID = elementIdx[axisMajor] + elementIdx[axisMedium] * numCartesianCells[axisMajor]
                + elementIdx[axisMinor] * numCartesianCells[axisMajor] * numCartesianCells[axisMedium];
            return firstElementID + numSimplicesPerCartesianCell * cartesianCellID + elementIdx[dim];
        }

        public int[] GetElementIdx(int elementID)
        {
            CheckElementID(elementID);

            var elementIdx = new int[dim + 1];
            elementIdx[dim] = (elementID - firstElementID) % numSimplicesPerCartesianCell;
            int cartesianCellID = (elementID - firstElementID) / numSimplicesPerCartesianCell;

            int numElementsPlane = numCartesianCells[axisMajor] * numCartesianCells[axisMedium];
            int mod = cartesianCellID % numElementsPlane;

            elementIdx[axisMinor] = cartesianCellID / numElementsPlane;
            elementIdx[axisMedium] = mod / numCartesianCells[axisMajor];
            elementIdx[axisMajor] = mod % numCartesianCells[axisMajor];
            return elementIdx;
        }

        public int[] GetElementConnectivity(int[] elementIdx)
        {
            CheckElementIdx(elementIdx);

            var nodeIDs = new int[numNodesPerElement];
            var nodeIdx = new int[dim]; // Avoid allocating an array per node
            for (int n = 0; n < numNodesPerElement; ++n)
            {
                int[][] offsetsOfTriangleNodes = elementNodeIdxOffsets[elementIdx[dim]];
                int[] offsetOfNode = offsetsOfTriangleNodes[n];
                nodeIdx[0] = elementIdx[0] + offsetOfNode[0];
                nodeIdx[1] = elementIdx[1] + offsetOfNode[1];
                nodeIdx[2] = elementIdx[2] + offsetOfNode[2];
                nodeIDs[n] = GetNodeID(nodeIdx);
            }

            return nodeIDs;
        }

        public int[] GetElementConnectivity(int elementID) => GetElementConnectivity(GetElementIdx(elementID));

        [Conditional("DEBUG")]
        private void CheckElementIdx(int[] elementIdx)
        {
            if (elementIdx.Length != dim + 1)
            {
                throw new ArgumentException($"Element index must be an array with length = {dim + 1}");
            }
            for (int d = 0; d < dim; ++d)
            {
                if ((elementIdx[d] < 0) || (elementIdx[d] >= numCartesianCells[d]))
                {
                    throw new ArgumentException($"Element index along dimension {d} must belong in [0, {numCartesianCells[d]})");
                }
            }
            if ((elementIdx[dim] < 0) || (elementIdx[dim] >= numSimplicesPerCartesianCell))
            {
                throw new ArgumentException($"Element index entry {dim} must belong in [0, {numSimplicesPerCartesianCell})");
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
                throw new ArgumentException($"Node index must be an array with length = {dim}");
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
            private readonly int[] numNodes;
            private int axisMajorChoice;
            private int axisMinorChoice;
            private int firstElementID;
            private int firstNodeID;
            private bool isAxisMajorMinorDefault;

            /// <summary>
            /// 
            /// </summary>
            /// <param name="minCoordinates"></param>
            /// <param name="maxCoordinates"></param>
            /// <param name="numNodes">Array with 3 positive integers.</param>
            public Builder(double[] minCoordinates, double[] maxCoordinates, int[] numNodes)
            {
                this.coordsMin = minCoordinates.Copy();
                this.coordsMax = maxCoordinates.Copy();
                this.numNodes = numNodes.Copy();

                // Defaults
                isAxisMajorMinorDefault = true;
                axisMajorChoice = int.MinValue;
                axisMinorChoice = int.MinValue;
                firstNodeID = 0;
                firstElementID = 0;
            }

            public UniformSimplicialMesh3D BuildMesh()
            {
                Validate();
                (int majorAxis, int mediumAxis, int minorAxis) = ChooseAxes();
                ValidateAxes(majorAxis, mediumAxis, minorAxis);

                return new UniformSimplicialMesh3D(coordsMin.Copy(), coordsMax.Copy(), numNodes.Copy(),
                    majorAxis, mediumAxis, minorAxis, firstNodeID, firstElementID);
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
                    entries.Add((numNodes[0], 0));
                    entries.Add((numNodes[1], 1));
                    entries.Add((numNodes[2], 2));
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
                if (numNodes.Length != dim) throw new ArgumentException($"Length of number of nodes must be {dim}.");
                for (int d = 0; d < dim; ++d)
                {
                    if (numNodes[d] < 1)
                    {
                        throw new ArgumentException($"Along axis {d}, there must be at least 1 node.");
                    }
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
