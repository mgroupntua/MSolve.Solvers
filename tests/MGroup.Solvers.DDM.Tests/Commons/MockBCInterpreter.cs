namespace MGroup.Solvers.DDM.Tests.Commons
{
	using System;
	using System.Collections.Generic;

	using MGroup.Constitutive.Structural;
	using MGroup.Constitutive.Structural.BoundaryConditions;
	using MGroup.Constitutive.Thermal;
	using MGroup.MSolve.Discretization.BoundaryConditions;
	using MGroup.MSolve.Discretization.Dofs;
	using MGroup.MSolve.Discretization.Entities;
	using MGroup.MSolve.Solution.AlgebraicModel;

	public class MockBCInterpreter : IAlgebraicModelInterpreter
	{
		private readonly IModel model;

		public MockBCInterpreter(IModel model)
		{
			ActiveDofs = new ActiveDofs();
			ActiveDofs.AddDof(StructuralDof.TranslationX);
			ActiveDofs.AddDof(StructuralDof.TranslationY);
			ActiveDofs.AddDof(StructuralDof.TranslationZ);
			ActiveDofs.AddDof(StructuralDof.RotationX);
			ActiveDofs.AddDof(StructuralDof.RotationY);
			ActiveDofs.AddDof(StructuralDof.RotationZ);
			ActiveDofs.AddDof(ThermalDof.Temperature);
			this.model = model;
		}

		public ActiveDofs ActiveDofs { get; }

		public IDictionary<(int, IDofType), (int, INode, double)> GetDirichletBoundaryConditionsWithNumbering() =>
			model.EnumerateSubdomains()
				.Select(x => new Tuple<int, IEnumerable<IBoundaryConditionSet<IDofType>>>(x.ID, model.EnumerateBoundaryConditions(x.ID)))
					.Select(i => i.Item2
						.SelectMany(x => x.EnumerateNodalBoundaryConditions(model.EnumerateElements(i.Item1)))
						.OfType<INodalDirichletBoundaryCondition<IDofType>>())
					.SelectMany(x => x)
					.OrderBy(x => x.Node.ID)
					.GroupBy(x => (x.Node.ID, x.DOF))
					.Select((x, Index) => (x.First().Node, (IDofType)x.Key.DOF, Index, x.Sum(a => a.Amount)))
				.ToDictionary(x => (x.Node.ID, x.Item2), x => (x.Index, x.Node, x.Item4));

		public IDictionary<(int, IDofType), (int, INode, double)> GetDirichletBoundaryConditionsWithNumbering(int subdomainID) =>
			model.EnumerateBoundaryConditions(subdomainID)
				.SelectMany(x => x.EnumerateNodalBoundaryConditions(model.EnumerateElements(subdomainID))).OfType<INodalDirichletBoundaryCondition<IDofType>>()
				.OrderBy(x => x.Node.ID)
				.GroupBy(x => (x.Node.ID, x.DOF))
				.Select((x, Index) => (x.First().Node, (IDofType)x.Key.DOF, Index, x.Sum(a => a.Amount)))
				.ToDictionary(x => (x.Node.ID, x.Item2), x => (x.Index, x.Node, x.Item4));
	}
}
