using System.Diagnostics;
using MGroup.MSolve.Discretization;

//TODO: Allow client to control which triangles to use and the order of their nodes.
//TODO: Element specific methods need testing
namespace MGroup.MSolve.Meshes.Structured
{
    /// <summary>
    /// Structured mesh with uniform distances between points on a cartesian grid. The elements are Tri3 (simplices) generated
    /// by dividing each Quad4 cell of an equivalent cartesian mesh into 2 subtriangles. The element indices are thus int[3], 
    /// where the first 2 entries are the index of the enclosing Quad4 of an equivalent cartesian mesh and the third entry is 
    /// 0 or 1 to denote one of the 2 subtriangles of said Quad4.
    /// </summary>
    public class UniformSimplicialMesh2D : IStructuredMesh
    {
        private const int dim = 2;
        private const int numNodesPerElement = 3;
        private const int numSimplicesPerCartesianCell = 2;

        private readonly int axisMajor;
        private readonly int axisMinor;
        private readonly double[] dx;
        private readonly int firstElementID;
        private readonly int firstNodeID;
        private readonly int[] numCartesianCells;

        /// <summary>
        /// Given the index (i,j) of the first node of a Quad4 with coordinates=(-1,-1), this list contains, for each subtriangle 
        /// of the Quad, the index offsets relative to (i,j) of the nodes of the subtriangle. The offset of node (i,j) is (0,0),
        /// the offset of (i+1, j) is (1,0) and so forth.
        /// </summary>
        private readonly List<int[][]> elementNodeIdxOffsets;

        private UniformSimplicialMesh2D(double[] minCoordinates, double[] maxCoordinates, int[] numNodes, int axisMajor,
            int axisMinor, int firstNodeID, int firstElementID)
        {
            this.MinCoordinates = minCoordinates;
            this.MaxCoordinates = maxCoordinates;
            this.NumNodes = numNodes;
            NumNodesTotal = NumNodes[0] * NumNodes[1];

            this.numCartesianCells = new int[dim];
            for (int d = 0; d < dim; ++d)
            {
                this.numCartesianCells[d] = numNodes[d] - 1;
            }
            NumElementsTotal = numSimplicesPerCartesianCell * numCartesianCells[0] * numCartesianCells[1];

            dx = new double[dim];
            for (int d = 0; d < dim; d++)
            {
                dx[d] = (maxCoordinates[d] - minCoordinates[d]) / numCartesianCells[d];
            }

            this.axisMajor = axisMajor;
            this.axisMinor = axisMinor;
            this.firstNodeID = firstNodeID;
            this.firstElementID = firstElementID;

            // The offsets of the nodes of the first triangle are
            elementNodeIdxOffsets = new List<int[][]>(numSimplicesPerCartesianCell);
            var triangle0NodeIdxOffsets = new int[numNodesPerElement][];
            triangle0NodeIdxOffsets[0] = new int[] { 0, 0 };
            triangle0NodeIdxOffsets[1] = new int[] { 1, 0 };
            triangle0NodeIdxOffsets[2] = new int[] { 0, 1 };
            elementNodeIdxOffsets.Add(triangle0NodeIdxOffsets);

            // And the offsets of the nodes of the second triangle are
            var triangle1NodeIdxOffsets = new int[numNodesPerElement][];
            triangle1NodeIdxOffsets[0] = new int[] { 1, 1 };
            triangle1NodeIdxOffsets[1] = new int[] { 0, 1 };
            triangle1NodeIdxOffsets[2] = new int[] { 1, 0 };
            elementNodeIdxOffsets.Add(triangle1NodeIdxOffsets);
        }

        public CellType CellType => CellType.Tri3;

        public int Dimension => dim;

        public double[] MinCoordinates { get; }

        public double[] MaxCoordinates { get; }

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
            for (int elementID = firstElementID; elementID < firstElementID + NumElementsTotal; ++elementID)
            {
                yield return (elementID, GetElementConnectivity(elementID));
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

        /// <summary>
        /// Example: Suppose a 3x4 mesh where there are 4 nodes with indices n7(1,2), n8(2,2), n10(1,3), n11(2,3) and 2 triangles 
        /// t10(n7,n8,n10), t11(n11,n10,n8). 
        /// t10 corresponds to <paramref name="elementIdx"/> = (1,2,0). 
        /// t11 corresponds to <paramref name="elementIdx"/> = (1,2,1). 
        /// </summary>
        /// <param name="elementIdx">The index of the element. An integer array with length = 3.</param>
        /// <returns></returns>
        public int GetElementID(int[] elementIdx)
        {
            CheckElementIdx(elementIdx);
            int cartesianCellID = elementIdx[axisMajor] + elementIdx[axisMinor] * numCartesianCells[axisMajor];
            return firstElementID + numSimplicesPerCartesianCell * cartesianCellID + elementIdx[dim];
        }

        public int[] GetElementIdx(int elementID) 
        {
            CheckElementID(elementID);
            var elementIdx = new int[dim + 1];
            elementIdx[dim] = (elementID - firstElementID) % numSimplicesPerCartesianCell;
            int cartesianCellID = (elementID - firstElementID) / numSimplicesPerCartesianCell; 
            elementIdx[axisMajor] = cartesianCellID % numCartesianCells[axisMajor];
            elementIdx[axisMinor] = cartesianCellID / numCartesianCells[axisMajor];
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
                    throw new ArgumentException($"Element index entry {d} must belong in [0, {numCartesianCells[d]})");
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
            private int firstElementID;
            private int firstNodeID;

            /// <summary>
            /// 
            /// </summary>
            /// <param name="minCoordinates"></param>
            /// <param name="maxCoordinates"></param>
            /// <param name="numElements">Array with 3 positive integers.</param>
            public Builder(double[] minCoordinates, double[] maxCoordinates, int[] numNodes)
            {
                this.coordsMin = minCoordinates.Copy();
                this.coordsMax = maxCoordinates.Copy();
                this.numNodes = numNodes.Copy();

                // Defaults
                axisMajorChoice = int.MinValue;
                firstNodeID = 0;
                firstElementID = 0;
            }

            public UniformSimplicialMesh2D BuildMesh()
            {
                Validate();
                (int majorAxis, int minorAxis) = ChooseAxes();
                return new UniformSimplicialMesh2D(coordsMin.Copy(), coordsMax.Copy(), numNodes.Copy(),
                    majorAxis, minorAxis, firstNodeID, firstElementID);
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
                    if (numNodes[0] <= numNodes[1])
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
                if (numNodes.Length != dim) throw new ArgumentException($"Length of number of nodes must be {dim}.");
                for (int d = 0; d < dim; ++d)
                {
                    if (numNodes[d] < 1)
                    {
                        throw new ArgumentException($"Along axis {d}, there must be at least 1 node.");
                    }
                }

                // Major axis
                if (axisMajorChoice != int.MinValue)
                {
                    if ((axisMajorChoice != 0) && (axisMajorChoice != 1)) throw new ArgumentException("Major axis must be 0 (x) or 1 (y)");
                }

                // First node & element IDs: no illegal values yet.
            }
        }
    }
}
