namespace MGroup.Solvers.DDM.DiscretizationExtensions
{
	using System;
	using System.Collections.Generic;
	using System.Linq;
	using System.Text;
	using System.Threading.Tasks;

	using MGroup.MSolve.Discretization;
	using MGroup.MSolve.Discretization.BoundaryConditions;
	using MGroup.MSolve.Discretization.Dofs;
	using MGroup.MSolve.Discretization.Entities;

	public static class ModelExtensions
	{
		public static INodalDirichletBoundaryCondition<IDofType>[] FindDirichletBCsOfNode(
			this IModel model, INode node, int subdomainID)
		{
			IEnumerable<IElementType> elements = model.EnumerateElements(subdomainID);
			var nodeConstraints = model.EnumerateBoundaryConditions(subdomainID)
					.Select(x => x.EnumerateNodalBoundaryConditions(elements))
					.OfType<INodalDirichletBoundaryCondition<IDofType>>()
					.Where(x => x.Node.ID == node.ID)
					.ToArray();
			return nodeConstraints;
		}

		public static IEnumerable<INodalDirichletBoundaryCondition<IDofType>> FindAllDirichletBCs(
			this IModel model, int subdomainID)
		{
			IEnumerable<IElementType> elements = model.EnumerateElements(subdomainID);
			var constraints = model.EnumerateBoundaryConditions(subdomainID)
				.SelectMany(x => x.EnumerateNodalBoundaryConditions(elements))
				.OfType<INodalDirichletBoundaryCondition<IDofType>>();
			return constraints;
		}
	}
}
