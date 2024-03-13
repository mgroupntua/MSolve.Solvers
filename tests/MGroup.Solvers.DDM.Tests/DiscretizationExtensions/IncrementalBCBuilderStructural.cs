namespace MGroup.Solvers.DDM.Tests.DiscretizationExtensions
{
	using System;
	using System.Collections.Generic;

	using MGroup.Constitutive.Structural;
	using MGroup.Constitutive.Structural.BoundaryConditions;
	using MGroup.MSolve.Discretization.Entities;
	using MGroup.Solvers.DDM.DiscretizationExtensions;

	//TODO: Similar functionality should be provided by the Model class itself.
	//TODO: Generic version that works independently of constitutive type
	//TODO: Make these work for distributed implementations of Model. So far they assume that if an entity exists, it is in memory.
	public class IncrementalBCBuilderStructural
	{
		private readonly TableExtension<INode, IStructuralDofType, double> dirichletBCs 
			= new TableExtension<INode, IStructuralDofType, double>();
		private readonly TableExtension<INode, IStructuralDofType, double> neumannBCs 
			= new TableExtension<INode, IStructuralDofType, double>();
		private readonly Model _model;
		private readonly bool throwExceptionIfDifferentBC;

		public IncrementalBCBuilderStructural(Model model)
		{
			_model = model;
		}

		/// <summary>
		/// Adds a new Dirichlet BC, only if there is no other BC at <paramref name="node"/>, <paramref name="dof"/>.
		/// </summary>
		public void AddDirichletBC(INode node, IStructuralDofType dof, double amount)
		{
			bool existsDirichlet = dirichletBCs.TryGetValue(node, dof, out double existingDirichlet);
			bool existsNeumann = neumannBCs.Contains(node, dof);
			if (existsDirichlet)
			{
				// Do nothing, unless the existing Dirichlet BC has a different value
				if (throwExceptionIfDifferentBC && (existingDirichlet != amount))
				{
					throw new Exception($"At node {node.ID}, dof = {dof}, u = {existingDirichlet} has already" +
						$" been prescribed. Cannot apply u = {amount} too.");
				}
			}
			else if (existsNeumann)
			{
				// Remove Neumann and keep the new Dirichlet
				neumannBCs.TryRemove(node, dof);
				dirichletBCs[node, dof] = amount;
			}
			else
			{
				dirichletBCs[node, dof] = amount;
			}
		}

		/// <summary>
		/// Adds a new Dirichlet BC, only if there is no other BC at <paramref name="node"/>, <paramref name="dof"/>.
		/// </summary>
		public void AddNeumannBC(INode node, IStructuralDofType dof, double amount)
		{
			bool existsDirichlet = dirichletBCs.Contains(node, dof);
			bool existsNeumann = neumannBCs.TryGetValue(node, dof, out double existingNeumann);
			if (existsDirichlet)
			{
				// Do nothing, Dirichlet overwrites any Neumann BC.
			}
			else if (existsNeumann)
			{
				// Sum the new and existing Neumann BCs.
				neumannBCs[node, dof] = amount + existingNeumann;
			}
			else
			{
				neumannBCs[node, dof] = amount;
			}
		}

		public void ConfirmBoundaryConditions()
		{
			var dirichlets = new List<INodalStructuralDirichletBoundaryCondition>();
			foreach ((INode node, IStructuralDofType dof, double amount) in dirichletBCs)
			{
				dirichlets.Add(new NodalDisplacement(node, dof, amount));
			}

			var neumanns = new List<INodalStructuralNeumannBoundaryCondition>();
			foreach ((INode node, IStructuralDofType dof, double amount) in neumannBCs)
			{
				neumanns.Add(new NodalLoad(node, dof, amount));
			}

			_model.BoundaryConditions.Add(new StructuralBoundaryConditionSet(dirichlets, neumanns));
		}
	}
}
