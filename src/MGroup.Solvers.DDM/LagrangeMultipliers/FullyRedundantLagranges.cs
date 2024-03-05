using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Text;

namespace MGroup.Solvers.DDM.LagrangeMultipliers
{
	public class FullyRedundantLagranges : ICrossPointStrategy
	{
		public List<(int subdomainPlus, int subdomainMinus)> ListSubdomainCombinations(IEnumerable<int> subdomainIDs)
		{
			var subdomains = new SortedSet<int>(subdomainIDs);
			int multiplicity = subdomains.Count;
			Debug.Assert(multiplicity > 1);
			int numCombos = (multiplicity * (multiplicity - 1)) / 2;
			var combos = new List<(int subdomainPlus, int subdomainMinus)>(numCombos);
			var processedSubdomains = new SortedSet<int>(subdomains);
			foreach (int subdomain1 in subdomains)
			{
				processedSubdomains.Remove(subdomain1);
				foreach (int subdomain2 in processedSubdomains)
				{
					// Lagrange multiplier between these 2 subdomains. The one with the min ID will be the positive one.
					if (subdomain1 < subdomain2)
					{
						combos.Add((subdomain1, subdomain2));
					}
					else
					{
						combos.Add((subdomain2, subdomain1));
					}
				}
			}
			return combos;
		}
	}
}
