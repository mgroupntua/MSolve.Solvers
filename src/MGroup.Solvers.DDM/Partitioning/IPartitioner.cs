namespace MGroup.Solvers.DDM.Partitioning
{
	using System.Collections.Generic;

	using MGroup.MSolve.Discretization.Entities;

	public interface IPartitioner
    {
        int NumSubdomainsTotal { get; }

        int GetClusterOfSubdomain(int subdomainID);

        IEnumerable<int> GetNeighboringSubdomains(int subdomainID);

        int GetSubdomainOfElement(int elementID);

        void Partition(IModel model);
    }
}
