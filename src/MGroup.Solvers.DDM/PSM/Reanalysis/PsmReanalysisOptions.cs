using System;
using System.Collections.Generic;
using System.Text;
using MGroup.Solvers.DDM.LinearSystem;

namespace MGroup.Solvers.DDM.PSM.Reanalysis
{
	public class PsmReanalysisOptions : ReanalysisOptions
	{
		public PsmReanalysisOptions(bool commonValueForFlags, IModifiedSubdomains modifiedSubdomains)
			: base(commonValueForFlags, modifiedSubdomains)
		{
			InterfaceProblemIndexer = commonValueForFlags;
			PreviousSolution = commonValueForFlags;
			RhsVectors = commonValueForFlags;
			Scaling = commonValueForFlags;
			SubdomainDofSubsets = commonValueForFlags;
			SubdomainSubmatrices = commonValueForFlags;
		}

		public bool InterfaceProblemIndexer { get; set; }

		public bool PreviousSolution { get; set; }

		public bool RhsVectors { get; set; }

		public bool Scaling { get; set; }

		public bool SubdomainDofSubsets { get; set; }

		public bool SubdomainSubmatrices { get; set; }

		public static new PsmReanalysisOptions CreateWithAllDisabled()
			=> new PsmReanalysisOptions(false, new NullModifiedSubdomains());

		public static new PsmReanalysisOptions CreateWithAllEnabled(IModifiedSubdomains modifiedSubdomains)
			=> new PsmReanalysisOptions(true, modifiedSubdomains);

	}
}
