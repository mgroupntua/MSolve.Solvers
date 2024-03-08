using MGroup.Solvers.DDM.LinearSystem;

namespace MGroup.Solvers.DDM.FetiDP
{
	public class FetiDPReanalysisOptions : ReanalysisOptions
	{
		public FetiDPReanalysisOptions(bool commonValueForFlags, IModifiedSubdomains modifiedSubdomains)
			: base(commonValueForFlags, modifiedSubdomains)
		{
			GlobalCoarseProblemDofs = commonValueForFlags;
			InterfaceProblemIndexer = commonValueForFlags;
			PreviousSolution = commonValueForFlags;
			RhsVectors = commonValueForFlags;
			Scaling = commonValueForFlags;
			SubdomainDofSubsets = commonValueForFlags;
			SubdomainSubmatrices = commonValueForFlags;
		}

		public bool GlobalCoarseProblemDofs { get; set; }

		public bool InterfaceProblemIndexer { get; set; }

		public bool PreviousSolution { get; set; }

		public bool RhsVectors { get; set; }

		public bool Scaling { get; set; }

		public bool SubdomainDofSubsets { get; set; }

		public bool SubdomainSubmatrices { get; set; }

		public static new FetiDPReanalysisOptions CreateWithAllDisabled()
			=> new FetiDPReanalysisOptions(false, new NullModifiedSubdomains());

		public static new FetiDPReanalysisOptions CreateWithAllEnabled(IModifiedSubdomains modifiedSubdomains)
			=> new FetiDPReanalysisOptions(true, modifiedSubdomains);
	}
}
