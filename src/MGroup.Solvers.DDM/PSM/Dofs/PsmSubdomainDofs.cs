using MGroup.MSolve.Discretization.Entities;
using MGroup.Solvers.DDM.Commons;
using MGroup.Solvers.DDM.LinearSystem;

//TODOMPI: DofTable should be replaced with an IntTable that stores ids, instead of actual references to nodes and dofs. 
//		This will make transfering it via MPI much faster.
//TODO: Naming convention for dofs (free/constrained, boundary/internal/corner/intercluster, subdomain/cluster/global) that will
//		be followed across all components
//TODOMPI: Replace DofTable with an equivalent class that uses integers. Also allow clients to choose sorted versions
namespace MGroup.Solvers.DDM.PSM.Dofs
{
	public class PsmSubdomainDofs
	{
		private readonly ISubdomain subdomain;
		private readonly ISubdomainLinearSystem linearSystem;

		//TODO: This is essential for testing and very useful for debugging, but not production code. Should I remove it?
		private readonly bool sortDofsWhenPossible;

		public PsmSubdomainDofs(ISubdomain subdomain, ISubdomainLinearSystem linearSystem, bool sortDofsWhenPossible = false)
		{
			this.subdomain = subdomain;
			this.linearSystem = linearSystem;
			this.sortDofsWhenPossible = sortDofsWhenPossible;
		}

		public IntDofTable DofOrderingBoundary { get; private set; }

		public int[] DofsBoundaryToFree { get; private set; }

		public int[] DofsInternalToFree { get; private set; }

		public bool IsEmpty => DofOrderingBoundary == null;

		public int NumFreeDofs { get; private set; }
		
		public void ReorderInternalDofs(DofPermutation permutation)
		{
			if (permutation.IsBetter)
			{
				DofsInternalToFree = permutation.ReorderKeysOfDofIndicesMap(DofsInternalToFree);
			}
		}

		/// <summary>
		/// Boundary/internal dofs
		/// </summary>
		public void SeparateFreeDofsIntoBoundaryAndInternal()
		{
			//TODOMPI: force sorting per node and dof
			var boundaryDofOrdering = new IntDofTable();
			var boundaryToFree = new List<int>();
			var internalToFree = new HashSet<int>();
			int subdomainBoundaryIdx = 0;

			IntDofTable freeDofs = linearSystem.DofOrdering.FreeDofs;
			IEnumerable<int> nodes = freeDofs.GetRows();
			if (sortDofsWhenPossible)
			{
				nodes = nodes.OrderBy(node => node);
			}

			foreach (int node in nodes) //TODO: Optimize access: Directly get INode, Dictionary<IDof, int>
			{
				IReadOnlyDictionary<int, int> dofsOfNode = freeDofs.GetDataOfRow(node);
				if (sortDofsWhenPossible)
				{
					var sortedDofsOfNode = new SortedDictionary<int, int>();
					foreach (var dofTypeIdxPair in dofsOfNode)
					{
						sortedDofsOfNode[dofTypeIdxPair.Key] = dofTypeIdxPair.Value;
					}
					dofsOfNode = sortedDofsOfNode;
				}

				if (subdomain.GetMultiplicityOfNode(node) > 1)
				{
					foreach (var dofTypeIdxPair in dofsOfNode)
					{
						int dofID = dofTypeIdxPair.Key;
						boundaryDofOrdering[node, dofID] = subdomainBoundaryIdx++;
						boundaryToFree.Add(dofTypeIdxPair.Value);
					}
				}
				else
				{
					foreach (var dofTypeIdxPair in dofsOfNode)
					{
						internalToFree.Add(dofTypeIdxPair.Value);
					}
				}
			}

			this.DofOrderingBoundary = boundaryDofOrdering;
			this.DofsBoundaryToFree = boundaryToFree.ToArray();
			this.DofsInternalToFree = internalToFree.ToArray();
		}
	}
}
