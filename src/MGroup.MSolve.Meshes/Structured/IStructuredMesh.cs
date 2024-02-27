using MGroup.MSolve.Discretization;

//TODO: Add non uniform rectilinear meshes, curvilinear (uniform or not) meshes, meshes with triangles/tetrahedra instead of 
//      quads/hexahedra, meshes with 2nd order and serendipity elements. Also rotated and skewed meshes. See this for ideas: 
//      https://axom.readthedocs.io/en/develop/axom/mint/docs/sphinx/sections/mesh_types.html#particlemesh
//TODO: Remove duplication between 2D/3D cartesian/simplicial and their builders.
namespace MGroup.MSolve.Meshes.Structured
{
    public interface IStructuredMesh
    {
        CellType CellType { get; }

        int Dimension { get; }

        double[] MinCoordinates { get; } //TODO: These probably apply only to unrotated & unskewed rectilinear meshes 

        double[] MaxCoordinates { get; }

        int NumElementsTotal { get; }

        int[] NumNodes { get; }

        int NumNodesPerElement { get; }

        int NumNodesTotal { get; }

        /// <summary>
        /// Enumerates the global IDs and coordinates of the nodes. The nodes are returned in increasing order of their id.
        /// </summary>
        IEnumerable<(int nodeID, double[] coordinates)> EnumerateNodes();

        /// <summary>
        /// For each element returns its global ID and the global IDs of its nodes. The elements are returned in increasing order
        /// of their id.
        /// </summary>
        IEnumerable<(int elementID, int[] nodeIDs)> EnumerateElements();

        int GetElementID(int[] elementIdx);

        int[] GetElementIdx(int elementID);

        int[] GetElementConnectivity(int[] elementIdx);

        int[] GetElementConnectivity(int elementID);

        int GetNodeID(int[] nodeIdx);

        int[] GetNodeIdx(int nodeID);

        double[] GetNodeCoordinates(int[] nodeIdx);

        double[] GetNodeCoordinates(int nodeID) => GetNodeCoordinates(GetNodeIdx(nodeID));

    }
}
