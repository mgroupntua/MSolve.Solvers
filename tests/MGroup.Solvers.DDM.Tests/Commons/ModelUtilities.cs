namespace MGroup.Solvers.DDM.Tests.Commons
{
	using MGroup.MSolve.Discretization;
	using MGroup.MSolve.Discretization.Entities;
	using MGroup.MSolve.Solution.AlgebraicModel;
	using MGroup.Solvers.DofOrdering;

	public static class ModelUtilities
	{
		public static ISubdomainFreeDofOrdering OrderDofs(ISubdomain subdomain, IAlgebraicModelInterpreter bcInterpreter)
		{
			var dofOrderer = new NodeMajorDofOrderingStrategy();
			(int numSubdomainFreeDofs, IntDofTable subdomainFreeDofs) = dofOrderer.OrderSubdomainDofs(subdomain, bcInterpreter);
			var dofOrdering = new SubdomainFreeDofOrderingCaching(
				numSubdomainFreeDofs, subdomainFreeDofs, bcInterpreter.ActiveDofs);
			return dofOrdering;
		}

		public static void DecomposeIntoSubdomains(this Model model, int numSubdomains, Func<int, int> getSubdomainOfElement)
		{
			model.SubdomainsDictionary.Clear();
			foreach (Node node in model.NodesDictionary.Values) node.Subdomains.Clear();
			foreach (IElementType element in model.ElementsDictionary.Values) element.SubdomainID = int.MinValue;

			for (int s = 0; s < numSubdomains; ++s)
			{
				model.SubdomainsDictionary[s] = new Subdomain(s);
			}
			foreach (IElementType element in model.ElementsDictionary.Values)
			{
				Subdomain subdomain = model.SubdomainsDictionary[getSubdomainOfElement(element.ID)];
				subdomain.Elements.Add(element);
			}

			model.ConnectDataStructures();
		}
	}
}
