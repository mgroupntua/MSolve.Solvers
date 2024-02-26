using System;
using System.Collections.Generic;
using System.Text;
using MGroup.LinearAlgebra.Matrices;
using MGroup.LinearAlgebra.Vectors;
using MGroup.Solvers.DDM.LinearAlgebraExtensions.Matrices;

//TODO: Use BLAS whenever possible
namespace MGroup.Solvers.DDM.LinearAlgebraExtensions
{
	public static class MatrixExtensions
	{
		/// <summary>
		/// <paramref name="matrix"/>[:, <paramref name="colIdx"/>] = <paramref name="matrix"/>[:,<paramref name="colIdx"/>] +
		/// <paramref name="wholeColumn"/>[:]
		/// </summary>
		/// <param name="matrix"></param>
		/// <param name="colIdx"></param>
		/// <param name="wholeColumn"></param>
		public static void AddColumn(this Matrix matrix, int colIdx, Vector wholeColumn)
			=> AxpyColumn(matrix, colIdx, +1.0, wholeColumn);

		/// <summary>
		/// <paramref name="matrix"/>[:, <paramref name="colIdx"/>] = <paramref name="matrix"/>[:,<paramref name="colIdx"/>] +
		/// <paramref name="colCoeff"/> * <paramref name="wholeColumn"/>[:]
		/// </summary>
		/// <param name="matrix"></param>
		/// <param name="colIdx"></param>
		/// <param name="colCoeff"></param>
		/// <param name="wholeColumn"></param>
		public static void AxpyColumn(this Matrix matrix, int colIdx, double colCoeff, Vector wholeColumn)
		{
			int start = colIdx * matrix.NumRows;
			double[] valuesMatrix = matrix.RawData;
			double[] valuesVector = wholeColumn.RawData;
			for (int i = 0; i < matrix.NumRows; ++i)
			{
				valuesMatrix[start + i] += colCoeff * valuesVector[i];
			}
		}

		public static void CopyFrom(this Matrix thisMatrix, Matrix otherMatrix)
		{
			Array.Copy(otherMatrix.RawData, thisMatrix.RawData, thisMatrix.RawData.Length);
		}

		/// <summary>
		/// <paramref name="matrix"/>[:, <paramref name="colIdx"/>] = <paramref name="matrix"/>[:,<paramref name="colIdx"/>] -
		/// <paramref name="wholeColumn"/>[:]
		/// </summary>
		/// <param name="matrix"></param>
		/// <param name="colIdx"></param>
		/// <param name="wholeColumn"></param>
		public static void SubtractColumn(this Matrix matrix, int colIdx, Vector wholeColumn)
			=> AxpyColumn(matrix, colIdx, -1.0, wholeColumn);
	}
}
