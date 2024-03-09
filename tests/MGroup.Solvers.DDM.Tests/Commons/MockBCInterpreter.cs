namespace MGroup.Solvers.DDM.Tests.Commons
{
	using System;
	using System.Collections.Generic;

	using MGroup.Constitutive.Structural;
	using MGroup.MSolve.Discretization.Dofs;
	using MGroup.MSolve.Discretization.Entities;
	using MGroup.MSolve.Solution.AlgebraicModel;

	public class MockBCInterpreter : IAlgebraicModelInterpreter
	{
		public MockBCInterpreter()
		{
			ActiveDofs = new ActiveDofs();
			ActiveDofs.AddDof(StructuralDof.TranslationX);
			ActiveDofs.AddDof(StructuralDof.TranslationY);
			ActiveDofs.AddDof(StructuralDof.TranslationZ);
			ActiveDofs.AddDof(StructuralDof.RotationX);
			ActiveDofs.AddDof(StructuralDof.RotationY);
			ActiveDofs.AddDof(StructuralDof.RotationZ);
		}

		public ActiveDofs ActiveDofs { get; }

		public IDictionary<(int NodeID, IDofType DOF), (int Index, INode Node, double Amount)> GetDirichletBoundaryConditionsWithNumbering() 
			=> throw new NotImplementedException();

		public IDictionary<(int NodeID, IDofType DOF), (int Index, INode Node, double Amount)> GetDirichletBoundaryConditionsWithNumbering(int subdomainID) 
			=> throw new NotImplementedException();
	}
}
