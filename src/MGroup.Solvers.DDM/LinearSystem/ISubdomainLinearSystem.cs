using System;
using System.Collections.Generic;
using System.Text;
using MGroup.LinearAlgebra.Matrices;
using MGroup.LinearAlgebra.Vectors;
using MGroup.Solvers.DofOrdering;

namespace MGroup.Solvers.DDM.LinearSystem
{
	public interface ISubdomainLinearSystem
	{
		ISubdomainFreeDofOrdering DofOrdering { get; }

		IMatrix Matrix { get; }

		Vector RhsVector { get; }

		Vector Solution { get; set; }

		int SubdomainID { get; }
	}
}
