namespace MGroup.Solvers.DDM.Partitioning
{
	using System.Diagnostics;

	using MGroup.MSolve.Discretization.Entities;
	using MGroup.MSolve.Discretization.Meshes.Structured;

	public class UniformMeshPartitioner2D : IPartitioner
    {
        private const int dim = 2;
        private const int invalidID = int.MinValue;

        private readonly UniformCartesianMesh2D mesh;
        private readonly int[] numClusters;
        private readonly int[] numSubdomains;
        private Dictionary<int, int> clustersOfSubdomains;
        private Dictionary<int, HashSet<int>> neighborsOfSubdomains;
        private Dictionary<int, int> subdomainsOfElements;

        public UniformMeshPartitioner2D(UniformCartesianMesh2D mesh, int[] numSubdomains, int[] numClusters)
        {
            this.mesh = mesh;
            this.numSubdomains = numSubdomains;
            this.numClusters = numClusters;
            this.NumSubdomainsTotal = numSubdomains[0] * numSubdomains[1];

            for (int d = 0; d < dim; ++d)
            {
                if (numSubdomains[d] % numClusters[d] != 0)
                {
                    throw new ArgumentException("The number of subdomains must be a multiple of the number of clusters per axis");
                }
                if (mesh.NumElements[d] % numSubdomains[d] != 0)
                {
                    throw new ArgumentException("The number of elements must be a multiple of the number of subdomains per axis");
                }
            }
        }

        public int NumSubdomainsTotal { get; }


        public int GetClusterOfSubdomain(int subdomainID) => clustersOfSubdomains[subdomainID];

        public IEnumerable<int> GetNeighboringSubdomains(int subdomainID) => neighborsOfSubdomains[subdomainID];

        public int GetSubdomainOfElement(int elementID) => subdomainsOfElements[elementID];

        public void Partition(IModel model)
        {
            PartitionElements();
            FindSubdomainNeighbors();
            ClusterSubdomains();
        }

        private void ClusterSubdomains()
        {
            var multiple = new int[dim];
            for (int d = 0; d < dim; ++d)
            {
                multiple[d] = numSubdomains[d] / numClusters[d];
            }

            clustersOfSubdomains = new Dictionary<int, int>();
            for (int sJ = 0; sJ < numSubdomains[1]; ++sJ)
            {
                int cJ = sJ / multiple[1];
                for (int sI = 0; sI < numSubdomains[0]; ++sI)
                {
                    int cI = sI / multiple[0];
                    int subdomainID = FindSubdomainID(new int[] { sI, sJ });
                    Debug.Assert(subdomainID != invalidID);
                    int clusterID = FindClusterID(new int[] { cI, cJ });
                    Debug.Assert(clusterID != invalidID);
                    clustersOfSubdomains[subdomainID] = clusterID;
                }
            }
        }

        private void FindSubdomainNeighbors()
        {
            neighborsOfSubdomains = new Dictionary<int, HashSet<int>>();
            for (int sJ = 0; sJ < numSubdomains[1]; ++sJ)
            {
                for (int sI = 0; sI < numSubdomains[0]; ++sI)
                {
                    int subdomainID = FindSubdomainID(new int[] { sI, sJ});
                    Debug.Assert(subdomainID != invalidID);
                    var neighborIDs = new HashSet<int>();

                    var neighborIndices = new int[8][];
                    neighborIndices[0] = new int[] { sI - 1, sJ - 1 }; 
                    neighborIndices[1] = new int[] { sI,     sJ - 1 };
                    neighborIndices[2] = new int[] { sI + 1, sJ - 1 };
                    neighborIndices[3] = new int[] { sI - 1, sJ     };
                    neighborIndices[4] = new int[] { sI + 1, sJ     };
                    neighborIndices[5] = new int[] { sI - 1, sJ + 1 };
                    neighborIndices[6] = new int[] { sI,     sJ + 1 };
                    neighborIndices[7] = new int[] { sI + 1, sJ + 1 };

                    foreach (int[] neighborIdx in neighborIndices)
                    {
                        int neighborID = FindSubdomainID(neighborIdx);

                        // For subdomains on the domain's boundary, some neighbors do not exist
                        if (neighborID != invalidID) 
                        {
                            neighborIDs.Add(neighborID);
                        }
                    }

                    neighborsOfSubdomains[subdomainID] = neighborIDs;
                }
            }
        }

        private void PartitionElements()
        {
            var multiple = new int[dim];
            for (int d = 0; d < dim; ++d)
            {
                multiple[d] = mesh.NumElements[d] / numSubdomains[d];
            }

            subdomainsOfElements = new Dictionary<int, int>();
            for (int eJ = 0; eJ < mesh.NumElements[1]; ++eJ)
            {
                int sJ = eJ / multiple[1];
                for (int eI = 0; eI < mesh.NumElements[0]; ++eI)
                {
                    int sI = eI / multiple[0];
                    int elementID = mesh.GetElementID(new int[] { eI, eJ });
                    int subdomainID = FindSubdomainID(new int[] { sI, sJ });
                    Debug.Assert(subdomainID != invalidID);
                    subdomainsOfElements[elementID] = subdomainID;
                }
            }
        }

        /// <summary>
        /// Returns <see cref="invalidID"/> if <paramref name="clusterIdx"/> is out of bounds.
        /// </summary>
        /// <param name="clusterIdx"></param>
        /// <returns></returns>
        private int FindClusterID(int[] clusterIdx)
        {
            for (int d = 0; d < dim; ++d)
            {
                if ((clusterIdx[d] < 0) || (clusterIdx[d] >= numClusters[d]))
                {
                    return invalidID;
                }
            }

            return clusterIdx[0] + clusterIdx[1] * numClusters[0];
        }

        /// <summary>
        /// Returns <see cref="invalidID"/> if <paramref name="subdomainIdx"/> is out of bounds.
        /// </summary>
        /// <param name="subdomainIdx"></param>
        /// <returns></returns>
        private int FindSubdomainID(int[] subdomainIdx)
        {
            for (int d = 0; d < dim; ++d)
            {
                if ((subdomainIdx[d] < 0) || (subdomainIdx[d] >= numSubdomains[d]))
                {
                    return invalidID;
                }
            }

            return subdomainIdx[0] + subdomainIdx[1] * numSubdomains[0];
        }
    }
}
