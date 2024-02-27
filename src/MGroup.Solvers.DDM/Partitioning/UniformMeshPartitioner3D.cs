using System.Diagnostics;
using MGroup.MSolve.Discretization.Entities;
using MGroup.MSolve.Meshes.Structured;

namespace MGroup.Solvers.DDM.Partitioning
{
	public class UniformMeshPartitioner3D : IPartitioner
	{
		private const int dim = 3;
		private const int invalidID = int.MinValue;

		private readonly UniformCartesianMesh3D mesh;
		private readonly int[] numClusters;
		private readonly int[][] numElementsPerSubdomainPerAxis;
		private readonly int[] numSubdomains;
		private readonly int[][] numSubdomainsPerClusterPerAxis;
		private Dictionary<int, int> clustersOfSubdomains;
		private Dictionary<int, HashSet<int>> neighborsOfSubdomains;
		private Dictionary<int, int> subdomainsOfElements;

		/// <summary>
		/// 
		/// </summary>
		/// <param name="mesh"></param>
		/// <param name="numSubdomains"></param>
		/// <param name="numClusters"></param>
		/// <param name="numElementsPerSubdomainPerAxis">
		/// For each axis d=0,1,2 <paramref name="numElementsPerSubdomainPerAxis"/> contains an int[] array. This array contains
		/// the number of elements of each subdomain. If this parameter is not provided: a) If the total number of elements along 
		/// an axis is a multiple of the number of subdomains per axis, then that multiplicity will be the number of elements 
		/// per subdomain along that axis. b) If it is not a multiple, then any extra elements will be distributed among a 
		/// subset of the subdomains, starting from the first subdomain of that axis.
		/// </param>
		public UniformMeshPartitioner3D(UniformCartesianMesh3D mesh, int[] numSubdomains, int[] numClusters, 
			int[][] numElementsPerSubdomainPerAxis = null)
		{
			this.mesh = mesh;
			this.numSubdomains = numSubdomains;
			this.numClusters = numClusters;
			this.NumSubdomainsTotal = numSubdomains[0] * numSubdomains[1] * numSubdomains[2];

			// Elements per subdomain
			if (numElementsPerSubdomainPerAxis != null)
			{
				// Check input
				for (int d = 0; d < dim; ++d)
				{
					int numElements = 0;
					if (numElementsPerSubdomainPerAxis[d].Length != numSubdomains[d])
					{
						throw new ArgumentException($"The number of subdomains along axis {d} is defined twice and differently.");
					}
					for (int s = 0; s < numElementsPerSubdomainPerAxis[d].Length; ++s)
					{
						if (numElementsPerSubdomainPerAxis[d][s] < 1)
						{
							throw new ArgumentException($"There must be at least 1 element per subdomain per axis.");
						}
						numElements += numElementsPerSubdomainPerAxis[d][s];
					}
					if (numElements != mesh.NumElements[d])
					{
						throw new ArgumentException($"The total number of elements in subdomains along axis {d} must be equal to" +
							$" the total number of elements in the mesh along that axis.");
					}
				}
				this.numElementsPerSubdomainPerAxis = numElementsPerSubdomainPerAxis;
			}
			else
			{
				this.numElementsPerSubdomainPerAxis = new int[dim][];
				for (int d = 0; d < dim; ++d)
				{
					this.numElementsPerSubdomainPerAxis[d] = new int[numSubdomains[d]];

					// All subdomains take at least numElements / numSubdomains elements.
					int div = mesh.NumElements[d] / numSubdomains[d];
					for (int s = 0; s < numSubdomains[d]; ++s)
					{
						this.numElementsPerSubdomainPerAxis[d][s] = div;
					}

					// Leftover elements are distributed among some subdomains
					int mod = mesh.NumElements[d] % numSubdomains[d];
					for (int s = 0; s < mod; ++s)
					{
						++this.numElementsPerSubdomainPerAxis[d][s];
					}
				}
			}

			// Subdomains per cluster
			this.numSubdomainsPerClusterPerAxis = new int[dim][];
			for (int d = 0; d < dim; ++d)
			{
				this.numSubdomainsPerClusterPerAxis[d] = new int[numClusters[d]];

				// All clusters take at least numSubdomains / numClusters subdomains.
				int div = numSubdomains[d] / numClusters[d];
				for (int c = 0; c < numClusters[d]; ++c)
				{
					this.numSubdomainsPerClusterPerAxis[d][c] = div;
				}

				// Leftover subdomains are distributed among some clusters (starting from the last ones)
				int mod = numSubdomains[d] % numClusters[d];
				for (int cc = 0; cc < mod; ++cc)
				{
					int c = numClusters[d] - 1 - cc;
					++this.numSubdomainsPerClusterPerAxis[d][c];
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
			clustersOfSubdomains = new Dictionary<int, int>();
			for (int sK = 0; sK < numSubdomains[2]; ++sK)
			{
				int cK = FindClusterOfSubdomain(2, sK);
				for (int sJ = 0; sJ < numSubdomains[1]; ++sJ)
				{
					int cJ = FindClusterOfSubdomain(1, sJ);
					for (int sI = 0; sI < numSubdomains[0]; ++sI)
					{
						int cI = FindClusterOfSubdomain(0, sI);
						int subdomainID = FindSubdomainID(new int[] { sI, sJ, sK });
						int clusterID = FindClusterID(new int[] { cI, cJ, cK });
						Debug.Assert(clusterID != invalidID);
						clustersOfSubdomains[subdomainID] = clusterID;
					}
				}
			}
		}

		private int FindClusterOfSubdomain(int axisIdx, int subdomainIdx)
		{
			int[] numSubdomainsPerCluster = numSubdomainsPerClusterPerAxis[axisIdx];
			int offset = 0;
			for (int c = 0; c < numSubdomainsPerCluster.Length; ++c)
			{
				if (subdomainIdx < offset + numSubdomainsPerCluster[c])
				{
					return c;
				}
				else
				{
					offset += numSubdomainsPerCluster[c];
				}
			}
			throw new ArgumentException($"Along axis {axisIdx} there are {numSubdomains[axisIdx]} subdomains." +
				$" Could not found requested subdomain with idx={subdomainIdx} along this axis");
		}

		private void FindSubdomainNeighbors()
		{
			neighborsOfSubdomains = new Dictionary<int, HashSet<int>>();
			for (int sK = 0; sK < numSubdomains[2]; ++sK)
			{
				for (int sJ = 0; sJ < numSubdomains[1]; ++sJ)
				{
					for (int sI = 0; sI < numSubdomains[0]; ++sI)
					{
						int subdomainID = FindSubdomainID(new int[] { sI, sJ, sK });
						Debug.Assert(subdomainID != invalidID);
						var neighborIDs = new HashSet<int>();

						var neighborIndices = new List<int[]>();
						for (int k = -1; k <= +1; ++k)
						{
							for (int j = -1; j <= +1; ++j)
							{
								for (int i = -1; i <= +1; ++i)
								{
									if ((i == 0) && (j == 0) && (k == 0))
									{
										continue; // Not a neighbor; the subdomain itself
									}
									neighborIndices.Add(new int[] { sI + i, sJ + j, sK + k });
								}
							}

						}

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
		}

		private int FindSubdomainOfElement(int axisIdx, int elementIdx)
		{
			int[] numElementsPerSubdomain = numElementsPerSubdomainPerAxis[axisIdx];
			int offset = 0;
			for (int s = 0; s < numElementsPerSubdomain.Length; ++s)
			{
				if (elementIdx < offset + numElementsPerSubdomain[s])
				{
					return s;
				}
				else
				{
					offset += numElementsPerSubdomain[s];
				}
			}
			throw new ArgumentException($"Along axis {axisIdx} there are {mesh.NumElements[axisIdx]} elements." +
				$" Could not found requested element with idx={elementIdx} along this axis");
		}

		private void PartitionElements()
		{
			subdomainsOfElements = new Dictionary<int, int>();
			for (int eK = 0; eK < mesh.NumElements[2]; ++eK)
			{
				int sK = FindSubdomainOfElement(2, eK);
				for (int eJ = 0; eJ < mesh.NumElements[1]; ++eJ)
				{
					int sJ = FindSubdomainOfElement(1, eJ);
					for (int eI = 0; eI < mesh.NumElements[0]; ++eI)
					{
						int sI = FindSubdomainOfElement(0, eI);
						int elementID = mesh.GetElementID(new int[] { eI, eJ, eK });
						int subdomainID = FindSubdomainID(new int[] { sI, sJ, sK });
						Debug.Assert(subdomainID != invalidID);
						subdomainsOfElements[elementID] = subdomainID;
					}
				}
			}
		}

		/// <summary>
		/// Returns <see cref="invalidID"/> if <paramref name="clusterIdx"/> is out of bounds.
		/// </summary>
		/// <param name="clusterIdx"></param>
		private int FindClusterID(int[] clusterIdx)
		{
			for (int d = 0; d < dim; ++d)
			{
				if ((clusterIdx[d] < 0) || (clusterIdx[d] >= numClusters[d]))
				{
					return invalidID;
				}
			}

			return clusterIdx[0] + clusterIdx[1] * numClusters[0] + clusterIdx[2] * numClusters[0] * numClusters[1];
		}

		/// <summary>
		/// Returns <see cref="invalidID"/> if <paramref name="subdomainIdx"/> is out of bounds.
		/// </summary>
		/// <param name="subdomainIdx"></param>
		private int FindSubdomainID(int[] subdomainIdx)
		{
			for (int d = 0; d < dim; ++d)
			{
				if ((subdomainIdx[d] < 0) || (subdomainIdx[d] >= numSubdomains[d]))
				{
					return invalidID;
				}
			}

			return subdomainIdx[0] + subdomainIdx[1] * numSubdomains[0] + subdomainIdx[2] * numSubdomains[0] * numSubdomains[1];
		}
	}
}
