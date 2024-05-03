namespace MGroup.Solvers.DDM
{
	using MGroup.Environments;
	using MGroup.LinearAlgebra.Distributed.Overlapping;
	using MGroup.MSolve.Discretization.Entities;
	using MGroup.Solvers.DofOrdering;

	public interface ISubdomainTopology
	{
		DistributedOverlappingIndexer CreateDistributedVectorIndexer(Func<int, IntDofTable> getSubdomainDofs);

		void FindCommonDofsBetweenSubdomains();

		void FindCommonNodesBetweenSubdomains();

		SortedSet<int> GetCommonNodesOfSubdomains(int localSubdomainID, int neighborSubdomainID);

		SortedSet<int> GetNeighborsOfSubdomain(int subdomainID);

		void Initialize(IComputeEnvironment environment, IModel model, Func<int, ISubdomainFreeDofOrdering> getSubdomainFreeDofs);

		DistributedOverlappingIndexer RecreateDistributedVectorIndexer(Func<int, IntDofTable> getSubdomainDofs,
			DistributedOverlappingIndexer previousIndexer, Func<int, bool> isModifiedSubdomain);

		void RefindCommonDofsBetweenSubdomains(Func<int, bool> isModifiedSubdomain);
	}
}
