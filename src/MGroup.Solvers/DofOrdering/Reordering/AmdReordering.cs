//TODO: The static factory methods are fine for now (and as defaults), but if the algorithm classes need options, then the user 
//      must instantiate them. There should probably be an IAmdAlgorithm interface. Or this class should not be named AMD.
namespace MGroup.Solvers.DofOrdering.Reordering
{
	using MGroup.LinearAlgebra.Implementations;
	using MGroup.LinearAlgebra.Reordering;
	using MGroup.MSolve.Discretization;
	using MGroup.MSolve.Discretization.Entities;

	/// <summary>
	/// Reorders the unconstrained freedom degrees according to the fill-reducing permutation calculated by the Approximate 
	/// Minimum Degree algorithm. Note that the pattern of the sparse matrix, i.e. the positions of its non-zero entries, must 
	/// be constructed and then passed to AMD. These might be costly operations and AMD might fail.
	/// </summary>
	public class AmdReordering : IDofReorderingStrategy
    {
        private readonly IReorderingAlgorithm amd;

		public AmdReordering(IImplementationProvider provider)
		{
			this.amd = new AmdSymmetricOrdering(provider);
		}

        public void ReorderDofs(ISubdomain subdomain, ISubdomainFreeDofOrdering originalOrdering)
        {
            originalOrdering.Reorder(amd, subdomain);
        }
    }
}
