using MGroup.LinearAlgebra.Distributed.IterativeMethods.PCG;
using MGroup.LinearAlgebra.Matrices;
using MGroup.Solvers.DDM.LinearSystem;

namespace MGroup.Solvers.DDM.FetiDP.InterfaceProblem
{
	public class ApproximatePcgResidualConvergence<TMatrix> : IPcgResidualConvergence
		where TMatrix : class, IMatrix
	{
		private readonly DistributedAlgebraicModel<TMatrix> algebraicModel;

		private double normF0;

		public ApproximatePcgResidualConvergence(DistributedAlgebraicModel<TMatrix> algebraicModel)
		{
			this.algebraicModel = algebraicModel;
		}

		public double EstimateResidualNormRatio(PcgAlgorithmBase pcg)
		{
			return pcg.PrecondResidual.Norm2() / normF0;
		}

		public void Initialize(PcgAlgorithmBase pcg)
		{
			normF0 = algebraicModel.LinearSystem.RhsVector.Norm2();
		}
	}
}
