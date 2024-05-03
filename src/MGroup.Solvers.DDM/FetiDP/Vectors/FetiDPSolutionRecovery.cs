namespace MGroup.Solvers.DDM.FetiDP.Vectors
{
	using System.Collections.Concurrent;

	using MGroup.Environments;
	using MGroup.LinearAlgebra.Distributed.Overlapping;
	using MGroup.LinearAlgebra.Matrices.Operators;
	using MGroup.LinearAlgebra.Vectors;
	using MGroup.Solvers.DDM.FetiDP.CoarseProblem;
	using MGroup.Solvers.DDM.FetiDP.Dofs;
	using MGroup.Solvers.DDM.FetiDP.Scaling;
	using MGroup.Solvers.DDM.FetiDP.StiffnessMatrices;

	public class FetiDPSolutionRecovery
	{
		private readonly IFetiDPCoarseProblem coarseProblem;
		private readonly IComputeEnvironment environment;

		private readonly Func<int, FetiDPSubdomainDofs> getSubdomainDofs;
		private readonly Func<int, SubdomainLagranges> getSubdomainLagranges;
		private readonly Func<int, IFetiDPSubdomainMatrixManager> getSubdomainMatrices;
		private readonly Func<int, FetiDPSubdomainRhsVectors> getSubdomainRhs;
		private readonly IFetiDPScaling scaling;

		public FetiDPSolutionRecovery(IComputeEnvironment environment, IFetiDPCoarseProblem coarseProblem, IFetiDPScaling scaling,
			Func<int, FetiDPSubdomainDofs> getSubdomainDofs, Func<int, SubdomainLagranges> getSubdomainLagranges,
			Func<int, IFetiDPSubdomainMatrixManager> getSubdomainMatrices, Func<int, FetiDPSubdomainRhsVectors> getSubdomainRhs)
		{
			this.environment = environment;
			this.coarseProblem = coarseProblem;
			this.scaling = scaling;
			this.getSubdomainDofs = getSubdomainDofs;
			this.getSubdomainLagranges = getSubdomainLagranges;
			this.getSubdomainMatrices = getSubdomainMatrices;
			this.getSubdomainRhs = getSubdomainRhs;
		}

		public void CalcPrimalSolution(DistributedOverlappingVector lagrangeSolution, 
			DistributedOverlappingVector primalSolutionAtFreeDofs)
		{
			// Calculate uc
			// yc[e] = fcCondensed[e] + Kcr[e] * inv(Krr[e]) * Dr[e]^T * lamda[e]
			// uc[e] = Ac * yc[e], Ac = Lc[e] * inv(KccCondensed[Global]) * Lc[e]^T
			var yce = new ConcurrentDictionary<int, Vector>();
			var uce = new ConcurrentDictionary<int, Vector>();
			var DreTimesLambda = new ConcurrentDictionary<int, Vector>();
			environment.DoPerNode(s =>
			{
				IFetiDPSubdomainMatrixManager subdomainMatrices = getSubdomainMatrices(s);
				SignedBooleanMatrixRowMajor Drs = getSubdomainLagranges(s).MatrixDr;
				Vector fcsCondensed = getSubdomainRhs(s).VectorFcCondensed;

				Vector DrsTimesLambda = Drs.Multiply(lagrangeSolution.LocalVectors[s], true);
				Vector ycs = subdomainMatrices.MultiplyKcrTimes(subdomainMatrices.MultiplyInverseKrrTimes(DrsTimesLambda));
				ycs.AddIntoThis(fcsCondensed);
				var xcs = Vector.CreateZero(ycs.Length);

				yce[s] = ycs;
				uce[s] = xcs;
				DreTimesLambda[s] = DrsTimesLambda;
			});
			coarseProblem.SolveCoarseProblem(yce, uce);

			// Calculate ur, uf
			environment.DoPerNode(s =>
			{
				// ur[s] = invKrr[e] * (fr[s] - Dr[s]^T * lamda[s] - Krc[s] * uc[s])
				IFetiDPSubdomainMatrixManager subdomainMatrices = getSubdomainMatrices(s);
				Vector ucs = uce[s];
				Vector temp = getSubdomainRhs(s).VectorFr.Copy();

				temp.SubtractIntoThis(DreTimesLambda[s]);
				temp.SubtractIntoThis(subdomainMatrices.MultiplyKrcTimes(ucs));
				Vector urs = subdomainMatrices.MultiplyInverseKrrTimes(temp);

				// Gather uc[s], ur[s] into uf[s]
				FetiDPSubdomainDofs subdomainDofs = getSubdomainDofs(s);
				var ufs = primalSolutionAtFreeDofs.LocalVectors[s];
				ufs.CopyNonContiguouslyFrom(subdomainDofs.DofsRemainderToFree, urs);
				ufs.CopyNonContiguouslyFrom(subdomainDofs.DofsCornerToFree, ucs);
			});

			// Average uf across all corresponding subdomains 
			environment.DoPerNode(s =>
			{
				var ufs = primalSolutionAtFreeDofs.LocalVectors[s];
				scaling.ScaleSubdomainRhsVector(s, ufs);
			});
			primalSolutionAtFreeDofs.SumOverlappingEntries();
		}
	}
}
