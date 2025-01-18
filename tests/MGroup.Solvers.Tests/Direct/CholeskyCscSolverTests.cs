namespace MGroup.Solvers.Tests.Direct
{
	using System;
	using System.Collections.Generic;
	using System.Linq;
	using System.Text;
	using System.Threading.Tasks;

	using MGroup.LinearAlgebra.Implementations;
	using MGroup.Solvers.Direct;
	using MGroup.Solvers.DofOrdering;
	using MGroup.Solvers.DofOrdering.Reordering;
	using MGroup.Solvers.Tests.Benchmarks;
	using Xunit;

	public class CholeskyCscSolverTests
	{
		[Theory]
		[MemberData(nameof(TestSettings.ProvidersToTest), MemberType = typeof(TestSettings))]
		public static void TestWithoutReordering(IImplementationProvider provider)
		{
			CantileverBeam benchmark = new CantileverBeam.Builder().BuildWithQuad4Elements(200, 10);

			var solverFactory = new CholeskyCscSolver.Factory(provider);
			solverFactory.DofOrderer = new DofOrderer(new NodeMajorDofOrderingStrategy(), new NullReordering());
			var algebraicModel = solverFactory.BuildAlgebraicModel(benchmark.Model);
			using (CholeskyCscSolver solver = solverFactory.BuildSolver(algebraicModel))
			{
				benchmark.RunAnalysisAndCheck(algebraicModel, solver);
			}
		}

		[Theory]
		[MemberData(nameof(TestSettings.ProvidersToTest), MemberType = typeof(TestSettings))]
		public static void TestWithAmdReordering(IImplementationProvider provider)
		{
			CantileverBeam benchmark = new CantileverBeam.Builder().BuildWithQuad4Elements(200, 10);

			var solverFactory = new CholeskyCscSolver.Factory(provider);
			//solverFactory.DofOrderer = new DofOrderer(new NodeMajorDofOrderingStrategy(), new AmdReordering(provider)); // default
			var algebraicModel = solverFactory.BuildAlgebraicModel(benchmark.Model);
			using (CholeskyCscSolver solver = solverFactory.BuildSolver(algebraicModel))
			{
				benchmark.RunAnalysisAndCheck(algebraicModel, solver);
			}
		}
	}
}
