//TODO: Needs a better name
namespace MGroup.Solvers.DDM
{
	using System.Diagnostics;

	/// <remarks>
	/// Similar to <see cref="SubdomainTopologyGeneral"/>, but assumes that all subdomains have the same dofs at their common 
	/// nodes, thus many operations can be done more efficiently and simply. 
	/// </remarks>
	public class SubdomainTopologyOptimized : SubdomainTopologyGeneral
	{
		public override void FindCommonDofsBetweenSubdomains()
		{
			environment.DoPerNode(subdomainID =>
			{
				Dictionary<int, DofSet> commonDofs = FindLocalSubdomainDofsAtCommonNodes(subdomainID);
				commonDofsBetweenSubdomains[subdomainID] = commonDofs;
			});
		}

		public override void RefindCommonDofsBetweenSubdomains(Func<int, bool> isModifiedSubdomain)
		{
			environment.DoPerNode(subdomainID =>
			{
				if (isModifiedSubdomain(subdomainID))
				{
					Dictionary<int, DofSet> commonDofs = FindLocalSubdomainDofsAtCommonNodes(subdomainID);
					commonDofsBetweenSubdomains[subdomainID] = commonDofs;
				}
				else
				{
					Debug.Assert(commonDofsBetweenSubdomains.ContainsKey(subdomainID));
				}
			});
		}
	}
}
