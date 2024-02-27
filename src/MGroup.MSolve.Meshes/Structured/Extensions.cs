using System;
using System.Collections.Generic;
using System.Text;

namespace MGroup.MSolve.Meshes.Structured
{
	public static class Extensions
	{
		public static bool AreContiguousUniqueIndices(this int[] numbers)
		{
			int[] copy = numbers.Copy();
			Array.Sort(copy);

			if (copy[0] != 0) return false;
			for (int i = 0; i < copy.Length - 1; ++i)
			{
				if (copy[i + 1] != copy[i] + 1) return false;
			}

			return true;
		}

		public static int[] Copy(this int[] array)
		{
			var copy = new int[array.Length];
			Array.Copy(array, copy, array.Length);
			return copy;
		}

		public static double[] Copy(this double[] array)
		{
			var copy = new double[array.Length];
			Array.Copy(array, copy, array.Length);
			return copy;
		}
	}
}
