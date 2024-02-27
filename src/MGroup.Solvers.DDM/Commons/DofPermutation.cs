using System;
using System.Collections.Generic;
using System.Text;

namespace MGroup.Solvers.DDM.Commons
{
	public class DofPermutation
	{
		private DofPermutation(bool isBetter, int[] permutationArray, bool permutationIsOldToNew)
		{
			this.IsBetter = isBetter;
			this.PermutationArray = permutationArray;
			this.PermutationIsOldToNew = permutationIsOldToNew;
		}

		public bool IsBetter { get; }
		public int[] PermutationArray { get; }
		public bool PermutationIsOldToNew { get; }

		public static DofPermutation Create(int[] permutationArray, bool permutationIsOldToNew)
			=> new DofPermutation(true, permutationArray, permutationIsOldToNew);

		public static DofPermutation CreateNoPermutation()
			=> new DofPermutation(false, null, false);

		public int[] ReorderKeysOfDofIndicesMap(int[] dofIndicesMap)
		{
			int numDofs = PermutationArray.Length;
			var result = new int[numDofs];
			if (PermutationIsOldToNew)
			{
				for (int i = 0; i < numDofs; ++i)
				{
					result[PermutationArray[i]] = dofIndicesMap[i]; // i is old index
				}
			}
			else
			{
				for (int i = 0; i < numDofs; ++i)
				{
					result[i] = dofIndicesMap[PermutationArray[i]]; // i is new index
				}
			}
			return result;
		}

		public int[] ReorderValuesOfDofIndicesMap(int[] dofIndicesMap)
		{
			int[] oldToNewValues;
			if (PermutationIsOldToNew)
			{
				oldToNewValues = PermutationArray;
			}
			else
			{
				oldToNewValues = new int[PermutationArray.Length];
				for (int i = 0; i < PermutationArray.Length; ++i)
				{
					oldToNewValues[PermutationArray[i]] = i;
				}
			}

			var result = new int[dofIndicesMap.Length];
			for (int i = 0; i < result.Length; ++i)
			{
				result[i] = oldToNewValues[dofIndicesMap[i]];
			}
			return result;
		}
	}
}
